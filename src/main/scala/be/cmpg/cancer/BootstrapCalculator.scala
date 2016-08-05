package be.cmpg.cancer

import scopt.OptionParser
import java.io.File
import scala.util.Random
import be.cmpg.graph.Network
import java.util.HashSet
import be.cmpg.graph.Gene
import scala.collection.mutable.HashMap
import be.cmpg.walk.SubNetworkSelector
import scala.collection.JavaConversions._
import be.cmpg.utils.StatUtils
import java.io.FileWriter
import org.json.JSONObject
import org.json.JSONArray
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import java.util.LinkedList
import java.io.PrintWriter
import scala.io.Source
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import java.io.FileOutputStream
import java.io.InputStreamReader

object BootstrapCalculator extends App {

  val helper = new CancerHelper

  val parser = helper.getBasicArgParser("SSA.ME Bootstraap calculation")
  parser.opt[Int]("bootstraapExperiments") action { (x, c) =>
    c.copy(bootstraapExperiments = x)
  } text ("Number of calculations with random inputs to calculate the bootstraap support. Large values increase processing time. (100 by default)")
  
  parser.opt[Double]("minBootstraapSupport") action { (x, c) =>
    c.copy(minBootstraapSupport = x)
  } text ("Minimum bootstrap support to be included in the network (0.9 by default)")
  
  parser.opt[Boolean]("useCGC") action { (x, c) =>
    c.copy(useCGC = x)
  } text ("Use the Census of Cancer Genes in the True Positives set (true by default)")

  parser.opt[Boolean]("useNCG") action { (x, c) =>
  c.copy(useNCG = x)
  } text ("Use the network of cancer genes in the True Positives set (true by default)")
  
  val configOpt = parser.parse(args, Config(

    /*
     * Default Values
     */
    iterations = 1000,
    reinforcement = 0.005,
    forgetfulness = 0.996))

  if (configOpt.isEmpty) {
    println("Error in parameters.")
    parser.showUsageAsError
    System.exit(1)
  }

  val config = configOpt.get

  val outputFileName = config.outputFolder+config.outputPrefix + "_bootstraap.randomization.tsv"
  val numberOfLines = {
    val outputFfile = new File(outputFileName)
    if (outputFfile.exists) Source.fromFile(outputFfile).getLines.size else 0
  }

  /**
   * Load results
   */
  val rankedGenes = new CSVReader(new FileReader(config.outputFolder+config.outputPrefix + "_nodes.tsv"), '\t')
    .readAll()
    .drop(1) // remove Header
    .map { fields => Gene(fields(0)) }.toList

  /**
   * For Coment 2 of Reviewer 1
   *
   * {
   * val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config)
   *
   * val writer = new PrintWriter(config.outputPrefix + "_bla.tsv")
   * genePatientMatrix.filter { case (key, p) => rankedGenes.contains(key.gene) }
   * .keys.groupBy { x => x.gene }
   * //.filter(_._2.size > 0)
   * .foreach { case (gene, values) => writer.println(gene.name + "\t" + values.size)}
   * writer.close()
   * println("bla "+config.outputPrefix + "_bla.tsv")
   * }
   *
   * Done.
   */
  val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config,false)

  if (numberOfLines < config.bootstraapExperiments) {
    /**
     * Load Original Data
     */


    val altered_genes = genePatientMatrix.map(_._1.gene).toSet
    
    /**
     * Load Gene Annotations
     */

    val genesToStudy = rankedGenes //(List("PIK3CA", "TP53", "MYC", "GAB2", "VAV2", "TTN", "UBC").map { x => Gene(x) } ++ rankedGenes)
    val writer = new PrintWriter(new FileOutputStream(new File(outputFileName), true))

    if (numberOfLines == 0) {
      /*
       * Create header if file is new
       */
      genesToStudy.foreach { gene => writer.print("\t" + gene.name) }
      writer.println()
      writer.flush()
    }

    val alterationsBySample = genePatientMatrix.toList.groupBy(_._1.sample)
    
    for (e <- numberOfLines to config.bootstraapExperiments) {

      val randomizedGenePatientMatrix = helper.bootstraapSamples(alterationsBySample)
      
      val boostraapNetworkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = randomizedGenePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical)

      val walkers = helper.buildWalkers(geneList.view.toSet, boostraapNetworkManager)

      println("*****************************")
      println("> rep: " + e + " ["+config.outputPrefix+"]")
      println("*****************************")
      
      boostraapNetworkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers))

      val bootstraapRankedGenes = boostraapNetworkManager.getRankedAllGenes().take(config.outputGenes).toList

      writer.print(e)
      genesToStudy.foreach { x => writer.print("\t" + bootstraapRankedGenes.indexOf(x)) }
      writer.println()
      writer.flush()
    }
    writer.close()

  }
  
  println("Loading randomization results...")

  val randomizationResults = new CSVReader(new FileReader(new File(outputFileName)), '\t').readAll()

  val statsByGene = {

    val statsByGene = randomizationResults.take(1).get(0).drop(1).map { gene => (Gene(gene), new LinkedList[Int]) }

    randomizationResults.drop(1).foreach { fields =>
      statsByGene.zip(fields.drop(1).map(_.toInt)).foreach { x => x._1._2.add(x._2) }
    }

    statsByGene.map( e => (e._1, e._2.view.toList)).toMap
  }
  
  /*rankedGenes.foreach { gene =>
    println(gene.name + "\t" + (statsByGene.getOrElse(gene, 0).toDouble / numberOfLines))
    writerbootstraap.println(gene.name + "\t" + (statsByGene.getOrElse(gene, 0).toDouble / numberOfLines))
  }*/
  
  
  println("Loading true positive gene sets...")
                                
  val malacardsGenes = if (config.otherGeneList==""){Set[Gene]()} else {Source.fromFile(config.otherGeneList).getLines.map(line => Gene(line.split("\t")(0))).toSet}
  
  val truePositiveSet =  (if (config.useCGC) helper.cgcGenes(config.inputFolder) else Set[Gene]()) ++ (if (config.useNCG) helper.ncgGenes(config.inputFolder) else Set[Gene]()) ++ malacardsGenes
  
  println("Calculating ppv with bootstraap support > "+ config.minBootstraapSupport +" ...")
                                      
  val writerbootstraap = new PrintWriter(config.outputFolder+config.outputPrefix + "_networksPPV.bootstraap.tsv")
  writerbootstraap.println("SubnetworkSize\tTruePositives\tFalsePositives\tPPV(%)\tgenes")
  
  val posibleSubnetworks = List.range(1, rankedGenes.size).map { subnetworkSize =>
    
    val genesToConsider = rankedGenes.take(subnetworkSize)
    
    val (_TP, _FP) = genesToConsider.map  { gene => 
                                            val geneStats = statsByGene.getOrElse(gene, List[Int]())
                                            val bootstraapSupport = geneStats.count { x => x >= 0 && x <= subnetworkSize }.toDouble / geneStats.size()
                                            if (bootstraapSupport >= config.minBootstraapSupport) Some(gene) else None
                                          }
                                          .flatten
                                          .partition { gene => truePositiveSet contains gene }
    
    writerbootstraap.println(subnetworkSize +"\t"+ _TP.size +"\t"+ _FP.size + "\t"+100*(_TP.size+1)/(_TP.size+_FP.size+1)+"\t"+(_TP++_FP).map { _.name }.mkString(","))
    
    (_TP, _FP, subnetworkSize)
  }
  writerbootstraap.close();

  
  val dummyNetworkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical)
  
  var selectedGenes = Set[Gene]()
  var rankCount = 0;
  //Map[Gene, Map[String, Any]]
  val additional = posibleSubnetworks.map { e =>
    if (selectedGenes.size < config.outputGenes) {
    	val genes = (e._1 ++ e._2).toSet
      val toReturn = (genes -- selectedGenes).map { gene => (gene, Map("rank" -> rankCount, "pvalue" -> 9))}
      
      if (!toReturn.isEmpty) {
        rankCount += 1
        selectedGenes ++= genes
        toReturn
      } else {
        Set.empty[(Gene, Map[String, Int])]
      }
    } else {
      Set.empty[(Gene, Map[String, Int])]
    }
  }.flatten.toMap
  
  MutualExclusivityPrintPattern.printPattern(config.outputPrefix+".selected", rankedGenes.filter { selectedGenes contains _ }, dummyNetworkManager, genePatientMatrix, config.outputFolder, config.inputFolder, additional,malacardsGenes)
  
  println("Output in " + config.outputFolder+config.outputPrefix + "_networksPPV.bootstraap.tsv")
  println("Network in " + config.outputFolder+config.outputPrefix + ".best_network.html")

  
  
  
}