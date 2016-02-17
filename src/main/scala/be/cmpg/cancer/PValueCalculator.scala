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

object PValueCalculator extends App {

  val random = new Random(System.nanoTime())

  val helper = new CancerHelper

  val parser = helper.getBasicArgParser("SSA.ME P-Value calculation")
  parser.opt[Int]("pvalueExperiments") action { (x, c) =>
    c.copy(pvalueExperiments = x)
  } text ("Number of calculations with random inputs to calculate the experimental p-values. Large values increase processing time. (1000 by default)")

  parser.opt[Int]("randomDistance") action { (x, c) =>
    c.copy(randomDistance = x)
  } text ("Distance to randomize the mutations. Each mutation will be assigned to a gene in its chromosomal neighborhood. (50000 bp by default)")

  val configOpt = parser.parse(args, Config(

    /*
     * Default Values
     */
    iterations = 1000,
    reinforcement = 0.005,
    forgetfulness = 0.996,
    refNetwork = Seq("HT", "hiII14", "reactome"),
    useRank = true,
    hypergeometricTest = false,
    minMutPerGene = 3,
    seedGenesMutations = 1,
    outputPrefix = "ME",
    pvalueExperiments = 1000,
    randomDistance = 50000))

  if (configOpt.isEmpty) {
    parser.showUsageAsError
    System.exit(1)
  }

  val config = configOpt.get

  val outputFfile = new File(config.outputPrefix + "_pvalue.randomization.tsv")

  val numberOfLines = if(outputFfile.exists) Source.fromFile(outputFfile).getLines.size else 0

  /**
   * Load results
   */
  val rankedGenes = new CSVReader(new FileReader(config.outputPrefix + "_nodes.tsv"), '\t')
    .readAll()
    .drop(1) // remove Header
    .map { fields => Gene(fields(0)) }.toList

  if (numberOfLines < config.pvalueExperiments) {
    /**
     * Load Original Data
     */

    val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config)

    /**
     * Load Gene Annotations
     */

    val geneAnnotations = new CSVReader(new FileReader("h19Genes.txt"), '\t')
      .readAll()
      .drop(1)
      .map { fields => GeneAnnotation(Gene(fields(0)), fields(1), fields(2).toInt, fields(3).toInt) }

    val annotationMap = geneAnnotations.map { ga =>
      val closeGenes = geneAnnotations.filter { other => ga.chrm == other.chrm && (other.end > ga.start - config.randomDistance) && (other.start < ga.end + config.randomDistance) }
      val sumGenesSize = closeGenes.foldLeft(0)((x, y) => x + y.size)

      (ga.gene, AnnotationInfo(closeGenes, sumGenesSize))
    }.toMap

    val genesToStudy = (List("PIK3CA", "TP53", "MYC", "GAB2", "VAV2", "TTN", "UBC").map { x => Gene(x) } ++ rankedGenes).toSet

    val writer = new PrintWriter(new FileOutputStream(outputFfile, true))

    if (numberOfLines == 0) {
      /*
       * Create header if file is new
       */
      genesToStudy.foreach { gene => writer.print("\t" + gene.name) }
      writer.println()
    }

    for (e <- numberOfLines to config.pvalueExperiments) {

      val pvalNetworkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = getRandomizedGenePatientMatrix(),
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        hypergeometricTest = config.hypergeometricTest)

      val walkers = helper.buildWalkers(geneList.view.toSet, pvalNetworkManager)
      
      println("*****************************")
      println("> rep: "+e)
      println("*****************************")

      pvalNetworkManager.run(config.iterations, geneList.view.toSet, Runtime.getRuntime().availableProcessors() / 2, defwalkers = Some(walkers))

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

    def getRandomizedGenePatientMatrix(): Map[PolimorphismKey, Polimorphism] = {

      val toReturn = new HashMap[PolimorphismKey, Polimorphism]

      genePatientMatrix.keys.map { key =>

        if (annotationMap.contains(key.gene)) {
          val annotationInfo = annotationMap(key.gene)
          var newGene: Gene = annotationInfo.closeGenes(0).gene
          var ran = random.nextInt(annotationInfo.sumGenesSize)

          annotationInfo.closeGenes
            .takeWhile(_ => ran > 0)
            .foreach { annotation =>
              ran -= annotation.size
              newGene = annotation.gene
            }

          toReturn.put(PolimorphismKey(newGene, key.sample), Polimorphism(newGene.name))
        }
      }

      toReturn.view.toMap
    }
  }

  val randomizationResults = new CSVReader(new FileReader(outputFfile), '\t').readAll()
  
  val statsByGene = {
    
    val geneStatistics = randomizationResults.take(1).get(0).drop(1).map { gene => (Gene(gene), new LinkedList[Int]) }

    randomizationResults.drop(1).foreach { fields =>
      //if (geneStatistics.size == fields.size + 1) {
        geneStatistics.zip(fields.drop(1).map(_.toInt)).foreach { x => x._1._2.add(x._2) }
      //} else {
        //throw new RuntimeException("ERROR: Number of genes in node file and randomization files are different.")
      //}
    }
  
    geneStatistics.map { case (gene, values) =>
      val indexInRealData = rankedGenes.indexOf(gene)
      (gene, values.count { x => x >= 0 && indexInRealData >= x}) // less than zero means not found
    }.toMap
  }
  val writer = new PrintWriter(config.outputPrefix + "_nodes.pvalues.tsv")
  rankedGenes.foreach { gene => 
    println(gene.name + "\t" + (statsByGene.getOrElse(gene, 0).toDouble / numberOfLines))
    writer.println(gene.name + "\t" + (statsByGene.getOrElse(gene, 0).toDouble / numberOfLines))
  }
  writer.close();
  
  class AnnotationInfo(geneAnn: GeneAnnotation) {
    val closeGenes = new LinkedList[GeneAnnotation]
  }

}

case class GeneAnnotation(gene: Gene, chrm: String, start: Int, end: Int) {
  lazy val size = end - start
}
case class AnnotationInfo(closeGenes: Seq[GeneAnnotation], sumGenesSize: Int)