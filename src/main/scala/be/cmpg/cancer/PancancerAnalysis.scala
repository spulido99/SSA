package be.cmpg.cancer

import be.cmpg.graph.NetworkReader
import java.nio.file.Paths
import be.cmpg.graph.Interaction
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.JavaConversions._
import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import scala.collection.mutable.HashSet
import be.cmpg.graph.Gene
import scala.collection.mutable.HashMap
import java.io.File
import be.cmpg.walk.SubNetworkSelector
import java.io.PrintWriter
import scala.collection.mutable.Buffer
import scala.collection.mutable.ArrayBuffer
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

object PancancerAnalysis extends App {

  case class Config(
    iterations: Int,
    reinforcement: Double,
    forgetfulness: Double,
    refNetwork: Seq[String],
    useRank: Boolean,
    outputPrefix: String,
    inputFolder: String = "",
    convergence: Double = 0.0,
    outputFolder: String = "",
    files: Seq[File] = Seq(),
    debug: Seq[String] = Seq())

  val parser = new scopt.OptionParser[Config]("SSA.ME.") {
    head("Mutual exclusivity analysis by Small Subnetwork Analysis", "0.1")
    opt[Int]('i', "iterations") action { (x, c) =>
      c.copy(iterations = x)
    } text ("iterations is the number of iterations (default: 5000)")

    opt[String]('b', "inputFolder") required() action { (x, c) =>
      c.copy(inputFolder = x)
    } text ("The folder which contains the networks and gene lists.")

    opt[String]('o', "outputFolder") action { (x, c) =>
      c.copy(outputFolder = x)
    } text ("The folder to where the output should be written.")

    opt[Double]('r', "reinforcement") action { (x, c) =>
      c.copy(reinforcement = x)
    } text ("The reinforncement for succesful learning (default: 0.0005)")

    opt[Double]('f', "forgetfulness") action { (x, c) =>
      c.copy(forgetfulness = x)
    } text ("The forgetfulness after each iteration (default: 0.9995)")

    opt[Double]('c', "convergence") action { (x, c) =>
      c.copy(convergence = x)
    } text ("The convergence required calculated as the RMSE between the posterior distributions (default: 0.0, i.e. finish by maxing iterations)")

    opt[String]('o', "outputPrefix") action { (x, c) =>
      c.copy(outputPrefix = x)
    } text ("The prefix for the analisis output files (default: ME)")

    opt[Seq[String]]('n', "refNetwork") action { (x, c) =>
      c.copy(refNetwork = x)
    } text ("The reference network (biogrid, kinase-substrate, encode or/and files. biogrid, kinase-substrate and encode used by default.)")

    /*    
    opt[Int]('g', "maxOutputGenes") action { (x, c) =>
      c.copy(maxOutputGenes = x)
    } text ("Jump by ")
    */

    opt[Boolean]('u', "useRank") action { (x, c) =>
      c.copy(useRank = x)
    } text ("Rank the gene scores instead of the scaled values (default: true)")

    opt[Seq[File]]('m', "files") required () action { (x, c) =>
      c.copy(files = x)
    } text ("Files with p-values for each gene")

    opt[Seq[String]]('d', "debug") action { (x, c) =>
      c.copy(debug = x)
    } text ("Debug mode. List of genes that want to be observed (e.g. -d TP53,MYC,")

    /*
    opt[Map[String, File]]('e', "expression") action { (x, c) =>
      c.copy(expression = x)
    } text ("expression files from GISTIC (... -e cnv_peaks=<file1>,exp=<file2>,cnv_thresholds=<file3> ...).")
    */
    help("help")

  }
  // parser.parse returns Option[C]

  val helper = new CancerHelper

  val configOpt = parser.parse(args, Config(

    /*
       * Default Values
       */
    iterations = 5000,
    reinforcement = 0.0005,
    forgetfulness = 0.9995,
    refNetwork = Seq("biogrid", "kinase-substrate", "encode"),
    useRank = false,
    outputPrefix = "PAN"))

  if (!configOpt.isDefined) {
    parser.showUsageAsError
    System.exit(1)
  }

  val config = configOpt.get

  val interactions = helper.loadNetwork(config.refNetwork, config.inputFolder)
  val network = new Network(interactions)
  val translateGenesToEntrez: Map[String, String] = network.genes.map(g => (g.name, g.name)).toMap

  var max = 0.0;
  val genePValues = new HashMap[Gene, PValueInfo]

  val origins = List("coding", "enhancer", "lncrna", "promoter")

  def parseOrigin(file: File): String = {
    for (origin <- origins) {
      if (file.getName.contains(origin))
        return origin
    }
    throw new RuntimeException
  }

  config.files.foreach { file =>

    val origin = parseOrigin(file)

    val enhancerfile = file.getName.contains("enhancer")
    val mutations = new CSVReader(new FileReader(file), '\t')
    var fields: Array[String] = mutations.readNext() // Drop first line
    do {
      fields = mutations.readNext()
      if (fields != null && fields.size > 5 && fields(5) != "NA") {

        val value = -Math.log(fields(5).toDouble)

        max = Math.max(max, value)

        if (!enhancerfile) {

          val gene = Gene(fields(0).split("::")(1))

          val pValueInfoOpt = genePValues.get(gene)
          if (pValueInfoOpt.isEmpty)
            genePValues += gene -> PValueInfo(value, origin)
          else if (pValueInfoOpt.get.max < value)
            genePValues += gene -> PValueInfo(value, origin)

        } else {
          fields(0)
            .split("::").drop(1)
            .map { _.split("_") }
            .filter { _.size > 1 }
            .map { x => Gene(x(1)) }
            .foreach { gene =>
              {
                val pValueInfoOpt = genePValues.get(gene)
                if (pValueInfoOpt.isEmpty)
                  genePValues += gene -> PValueInfo(value, origin)
                else if (pValueInfoOpt.get.max < value)
                  genePValues += gene -> PValueInfo(value, origin)

                if (gene.name == "ADAM12")
                  println(genePValues(gene))
              }
            }
        }
      }

    } while (fields != null)
  }

  println("Max Value: " + max)

  val networkManager = new PancancerNetworkManager(
    network = network,
    genePValues = genePValues.view.toMap,
    pheromone = config.reinforcement,
    evaporation = config.forgetfulness,
    ranked = config.useRank,

    convergenceThreshold = config.convergence)

  //val stats = new DescriptiveStatistics
  //genePValues.values.map { _.getMax }.foreach { stats.addValue(_) }
  //val threshold = stats.getPercentile(0.5)
  //val geneList = genePValues.filter { e => network.genes.contains(e._1) && e._2.getMax > threshold }.keySet
  val geneList = genePValues.filter { e => network.genes.contains(e._1) }.keySet

  println("Genes: " + genePValues.size)
  println("Genes(InNetwork): " + genePValues.filter { e => network.genes.contains(e._1) }.size)
  //println("Walker Threshold: "+threshold)
  println("Walkers: " + geneList.size)

  if (!config.debug.isEmpty)
    networkManager.debug = Some(config.debug.map { Gene(_) }.toSet, 20)
  networkManager.setUserFakeMax(max)

  networkManager.run(config.iterations, geneList, Runtime.getRuntime().availableProcessors() / 2, storeNodePHistory = true)

  //val rankedGenes = networkManager.rankedGenes

  /*
   * For printing results purposes:
   */
  /*
   * Real PANCANCER
   */
  /*val rankedGenes = List("TP53", "PAX5", "MYC", "CREBBP", "VHL", "ACTB", 
                          "CTNNB1", "NFATC2", "COL4A2", "COL4A1", "CDKN2A", 
                          "EP300", "DDX5", "PADI4", "BZRAP1", "NHP2L1", 
                          "SUMO2", "AKT1", "BRAF", "CDKN1A", "MDM2", "EGR1", 
                          "E2F1", "TAF1", "NFE2L2", "PRKACA", "SUZ12").map(Gene(_))*/

  /*
   * Random PANCACNER
   */
  //val rankedGenes = List("ADAM12", "IGFBP3").map(Gene(_))

  /*
   * Real BRCA
   */
  val rankedGenes = List("TP53", "NOC2L", "SUMO2", "PTEN", "AKT1", "MAP2K4", "PIK3CA", "PIK3R1", "GATA3", "AXIN1", "KRT18", "DVL1").map(Gene(_))

  val genesToReport = rankedGenes

  println("Genes printed in network: " + rankedGenes.size)

  println("Gene Rank: ")
  networkManager.getRankedAllGenes().take(15).foreach { gene =>
    {
      val node = networkManager.getNetwork().getNode(gene)
      println(gene.name + "\t" + genePValues.getOrElse(gene, PValueInfo(0.0, null)).max + "\t" + node.posteriorProbability)
    }
  }
  println("...")
  networkManager.getRankedAllGenes().takeRight(15).foreach { gene =>
    {
      val node = networkManager.getNetwork().getNode(gene)
      println(gene.name + "\t" + genePValues.getOrElse(gene, PValueInfo(0.0, null)).max + "\t" + node.posteriorProbability)
    }
  }

  val additionalNodeInfo = rankedGenes.map { g =>
    val info = genePValues.get(g)
    if (info.isDefined)
      g -> Map[String, Any]("pvalue" -> info.get.max, "origin" -> info.get.origin)
    else
      g -> Map[String, Any]("pvalue" -> 0.0, "origin" -> "unkown")
  }.toMap

  MutualExclusivityPrintPattern.printPattern(config.outputPrefix, rankedGenes.view.toList, networkManager, Map(), config.outputFolder, config.inputFolder, additionalNodeInfo)

  val historicPWriter = new PrintWriter("genesHistoricPs_" + config.reinforcement + "_" + config.forgetfulness + ".out")
  rankedGenes.foreach(g => {
    historicPWriter.print(g.name)
    val node = network.getNode(g)
    node.probabilityHistory.foreach(p => historicPWriter.print("\t" + (p * 100).toInt))
    historicPWriter.println()
  })
  historicPWriter.close()

  val modulesPWriter = new PrintWriter(config.outputFolder + config.outputPrefix + "_modules.txt")

  modulesPWriter.println("p-Value\tModule")
  networkManager.getNetwork().getNodes()
    .map(e => e.bestSubnetwork)
    .foreach { e =>
      println(math.exp(-e._2) + "\t" + e._1)
      modulesPWriter.println(math.exp(-e._2) + "\t" + e._1)
    }
  modulesPWriter.close()

}
case class PValueInfo(max: Double, origin: String)