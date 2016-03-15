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

object PValueCalculator extends App {

  val helper = new CancerHelper

  val parser = helper.getBasicArgParser("SSA.ME P-Value calculation")
  parser.opt[Int]("pvalueExperiments") action { (x, c) =>
    c.copy(pvalueExperiments = x)
  } text ("Number of calculations with random inputs to calculate the experimental p-values. Large values increase processing time. (1000 by default)")

  parser.opt[String]("pvaluetype") action { (x, c) =>
    c.copy(pvaluetype = x)
  } text ("P-value calculation to run. [distance, genenames, closegenes, bootstraap] (default: distance) ")

  parser.opt[Int]("closestGenes") action { (x, c) =>
    c.copy(closestGenes = x)
  } text ("[Used only for pvaluetype:closegenes] Number of genes to randomize the mutations. Each mutation will be assigned to a gene in its chromosomal neighborhood. (5 genes by default)")
  
  parser.opt[Int]("genomicdistnace") action { (x, c) =>
  c.copy(genomicdistnace = x)
  } text ("[Used only for pvaluetype:distance] Genomic distance to find genes to randomize the mutations. Each mutation will be assigned to a gene in its chromosomal neighborhood. (50000 bp by default)")
  

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

  val outputFileName = config.outputPrefix + "_pvalue."+config.pvaluetype+".randomization.tsv"
  val numberOfLines = {
    val outputFfile = new File(outputFileName)
    if (outputFfile.exists) Source.fromFile(outputFfile).getLines.size else 0
  }

  /**
   * Load results
   */
  val rankedGenes = new CSVReader(new FileReader(config.outputPrefix + "_nodes.tsv"), '\t')
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

  if (numberOfLines < config.pvalueExperiments) {
    /**
     * Load Original Data
     */

    val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config)

    val altered_genes = genePatientMatrix.map(_._1.gene).toSet
    
    /**
     * Load Gene Annotations
     */

    val annotationMap = if (config.pvaluetype == "distance" || config.pvaluetype == "closegenes") {
      
      println("Loading human gene annotation...")
      val geneAnnotations = new CSVReader(new InputStreamReader(getClass.getResourceAsStream("/h19Genes.txt")), '\t')
        .readAll()
        .drop(1)
        .map { fields => GeneAnnotation(Gene(fields(0)), fields(1), fields(2).toInt, fields(3).toInt) }
        .filter { ga => geneList contains ga.gene }
  
      geneAnnotations.map { ga =>
  
        val closeGenes = if (config.pvaluetype == "distance")
              geneAnnotations.filter { other => ga.chrm == other.chrm }.sortBy { other => math.min(math.abs(ga.start - other.end), math.abs(ga.end - other.start)) }.take(config.closestGenes)
            else if (config.pvaluetype == "closegenes")
              geneAnnotations.filter { other => ga.chrm == other.chrm && (other.end > ga.start - config.genomicdistnace) && (other.start < ga.end + config.genomicdistnace) }
            else
              throw new RuntimeException("p-value type not known: "+config.pvaluetype)
        val sumGenesSize = closeGenes.foldLeft(0)((x, y) => x + y.size)
  
        (ga.gene, AnnotationInfo(closeGenes, sumGenesSize))
      }.toMap
    } else 
      Map[Gene, AnnotationInfo]()

    val genesToStudy = rankedGenes //(List("PIK3CA", "TP53", "MYC", "GAB2", "VAV2", "TTN", "UBC").map { x => Gene(x) } ++ rankedGenes)

    val writer = new PrintWriter(new FileOutputStream(outputFileName, true))

    if (numberOfLines == 0) {
      /*
       * Create header if file is new
       */
      genesToStudy.foreach { gene => writer.print("\t" + gene.name) }
      writer.println()
    }

    val alterationsBySample = if (config.pvaluetype == "bootstraap") genePatientMatrix.toList.groupBy(_._1.sample) else Map[String, List[(PolimorphismKey, Polimorphism)]]()
    
    for (e <- numberOfLines to config.pvalueExperiments) {

      val randomizedGenePatientMatrix = if (config.pvaluetype == "distance" || config.pvaluetype == "closegenes") {
          helper.getCloseGenesRandomizedGenePatientMatrix(genePatientMatrix, annotationMap)
        } else if (config.pvaluetype == "bipartite") {
      	  helper.randomizeBipartiteGenePatientMatrix(genePatientMatrix, all_samples, altered_genes intersect network.genes)
        } else if (config.pvaluetype == "genenames") {
          helper.randomizeGeneNames(genePatientMatrix, altered_genes intersect network.genes)
        } else if (config.pvaluetype == "bootstraap") {
          helper.bootstraapSamples(alterationsBySample)
        } else {
          throw new RuntimeException("p-value type not known: "+config.pvaluetype)
        }
      
      val pvalNetworkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = randomizedGenePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical)

      val walkers = helper.buildWalkers(geneList.view.toSet, pvalNetworkManager)

      println("*****************************")
      println("> rep: " + e)
      println("*****************************")

      pvalNetworkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers))

      val pvalRankedGenes = pvalNetworkManager.getRankedAllGenes().take(config.outputGenes).toList

      writer.print(e)
      genesToStudy.foreach { x => writer.print("\t" + pvalRankedGenes.indexOf(x)) }
      writer.println()
      writer.flush()
      /*
      rankedGenes.zipWithIndex.foreach {
        case (gene: Gene, index: Int) =>
        //pvalRankedGenes.indexOf(elem)
      }
      *
      */

    }
    writer.close()

  }
  
  println("Loading randomization results...")

  val randomizationResults = new CSVReader(new FileReader(outputFileName), '\t').readAll()

  val statsByGene = {

    val geneStatistics = randomizationResults.take(1).get(0).drop(1).map { gene => (Gene(gene), new LinkedList[Int]) }

    randomizationResults.drop(1).foreach { fields =>
      geneStatistics.zip(fields.drop(1).map(_.toInt)).foreach { x => x._1._2.add(x._2) }
    }

    geneStatistics.map {
      case (gene, values) =>
        val indexInRealData = rankedGenes.indexOf(gene)
        (gene, values.count { x => x >= 0 && indexInRealData >= x }) // less than zero means not found
    }.toMap
  }
  val writerPvalues = new PrintWriter(config.outputPrefix + "_nodes.pvalues.tsv")
  rankedGenes.foreach { gene =>
    println(gene.name + "\t" + (statsByGene.getOrElse(gene, 0).toDouble / numberOfLines))
    writerPvalues.println(gene.name + "\t" + (statsByGene.getOrElse(gene, 0).toDouble / numberOfLines))
  }
  writerPvalues.close();

}

