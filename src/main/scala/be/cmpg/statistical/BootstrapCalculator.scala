package be.cmpg.statistical

import java.io.File
import be.cmpg.graph.Gene
import scala.collection.JavaConversions._
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import java.util.LinkedList
import java.io.PrintWriter
import scala.io.Source
import java.io.FileOutputStream
import be.cmpg.cancer.Config
import be.cmpg.cancer.ME.MutualExclusivityNetworkManager
import be.cmpg.cancer.PolymorphismKey
import be.cmpg.cancer.Polymorphism
import be.cmpg.cancer.ArgumentsParser
import be.cmpg.cancer.DataLoader
import be.cmpg.walk.WalkerFactory
import be.cmpg.cancer.PrintFiles

object BootstrapCalculator extends App {

  val parser = ArgumentsParser.getBasicArgParser("SSA.ME Bootstrap calculation")
  parser.opt[Int]("bootstrapExperiments") action { (x, c) =>
    c.copy(bootstrapExperiments = x)
  } text ("Number of calculations with random inputs to calculate the bootstrap support. Large values increase processing time. (100 by default)")

  parser.opt[Int]("bootstrapOutPutGenes") action { (x, c) =>
    c.copy(bootstrapOutPutGenes = x)
  } text ("Maximum number of genes to be shown in the output. Will be less if less genes get selected in total. 100 by default")

  parser.opt[Double]("minBootstrapSupport") action { (x, c) =>
    c.copy(minBootstrapSupport = x)
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

  val outputFileName = config.outputFolder + config.outputPrefix + "_bootstrap.randomization.tsv"
  val numberOfLines = {
    val outputFfile = new File(outputFileName)
    if (outputFfile.exists) Source.fromFile(outputFfile).getLines.size else 0
  }

  /**
   * Load results
   */
  val rankedGenes = new CSVReader(new FileReader(config.outputFolder + config.outputPrefix + "_nodes.tsv"), '\t')
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
  val (interactions, network, genePatientMatrix, all_samples, geneList) = DataLoader.loadData(config, false)

  if (numberOfLines < config.bootstrapExperiments) {
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

    for (e <- numberOfLines to config.bootstrapExperiments) {

      val randomizedGenePatientMatrix = DataRandomizations.bootstrapSamples(alterationsBySample)

      val bootstrapNetworkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = randomizedGenePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical)

      val walkers = WalkerFactory.buildHi2iiFungusWalker(geneList.view.toSet, bootstrapNetworkManager)

      println("*****************************")
      println("> rep: " + e + " [" + config.outputPrefix + "]")
      println("*****************************")

      bootstrapNetworkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers))

      val bootstrapRankedGenes = bootstrapNetworkManager.getRankedAllGenes().take(config.outputGenes).toList

      writer.print(e)
      genesToStudy.foreach { x => writer.print("\t" + bootstrapRankedGenes.indexOf(x)) }
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

    statsByGene.map(e => (e._1, e._2.view.toList)).toMap
  }

  println("Loading true positive gene sets...")

  val malacardsGenes = if (config.otherGeneList == "") { Set[Gene]() } else { Source.fromFile(config.otherGeneList).getLines.map(line => Gene(line.split("\t")(0))).toSet }

  val truePositiveSet = (if (config.useCGC) DataLoader.cgcGenes(config.inputFolder) else Set[Gene]()) ++ (if (config.useNCG) DataLoader.ncgGenes(config.inputFolder) else Set[Gene]()) ++ malacardsGenes

  println("Calculating ppv with bootstrap support > " + config.minBootstrapSupport + " ...")

  val writerbootstrap = new PrintWriter(config.outputFolder + config.outputPrefix + "_networksPPV.bootstrap.tsv")
  writerbootstrap.println("SubnetworkSize\tTruePositives\tFalsePositives\tPPV(%)\tgenes")

  // This "subnetworkSize" is not the size of the subnetwork as it is simply counting and multiple genes can be added or no genes can be added when going one step further. In order to determine the subnetwork size correctly following code should be used:

  // Should also add code to actually calculate 100 best genes instead of randomly from the previous 100 best. So the Node.tsv files should be enlarged and the first 100 should be considered here.

  var posibleSubnetworks = List[(List[Gene], List[Gene], Int)]()

  //  val posibleSubnetworks = List.range(1, config.bootstrapOutPutGenes).map { counter =>

  var stopCriterium = false
  var subnetworksWithDataToPrint = List[String]()

  var i = 0

  while (!stopCriterium) {

    //    val genesToConsider = try { rankedGenes.take(counter) } catch { case e: Exception => rankedGenes.take(rankedGenes.size) }

    if (i > rankedGenes.size) { stopCriterium = true }

    val genesToConsider = try { rankedGenes.take(i) } catch {
      case e: Exception => {
        rankedGenes.take(rankedGenes.size)
      }
    }

    val (_TP, _FP) = genesToConsider.map { gene =>
      val geneStats = statsByGene.getOrElse(gene, List[Int]())
      // Gene has to be selected in a specific bootstrap run (x >= 0) and has to be lower then the counter in order to only take 
      //      val bootstraapSupport = geneStats.count { x => x >= 0 && x <= counter }.toDouble / geneStats.size()
      val bootstrapSupport = geneStats.count { x => x >= 0 && x <= i }.toDouble / geneStats.size()
      if (bootstrapSupport >= config.minBootstrapSupport) Some(gene) else None
    }
      .flatten
      .partition { gene => truePositiveSet contains gene }

    val subnetworkSize = _TP.size + _FP.size
    val PPV = try { (100 * (_TP.size) / (_TP.size + _FP.size)) } catch { case e: Exception => 0 }

    subnetworksWithDataToPrint = subnetworksWithDataToPrint.:+(subnetworkSize + "\t" + _TP.size + "\t" + _FP.size + "\t" + PPV + "\t" + (_TP ++ _FP).map { _.name }.mkString(","))

    // The +1 in the PPV formula could indeed be used to guard against 0 in the denominator but it gives significant errors for low subnetwork size and a denominator of 0 can not occur in case of subnetworks (there is always at least one gene in the network which is either FP or TP)                                       
    //    writerbootstraap.println(subnetworkSize + "\t" + _TP.size + "\t" + _FP.size + "\t" + PPV + "\t" + (_TP ++ _FP).map { _.name }.mkString(","))

    posibleSubnetworks = posibleSubnetworks.:+(_TP, _FP, subnetworkSize)
    i = i + 1

    if (subnetworkSize > config.bootstrapOutPutGenes) { stopCriterium = true }
  }
  subnetworksWithDataToPrint.distinct.foreach(line => writerbootstrap.println(line))
  writerbootstrap.close();

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
    if (selectedGenes.size < config.bootstrapOutPutGenes) {
      val genes = (e._1 ++ e._2).toSet
      val toReturn = (genes -- selectedGenes).map { gene => (gene, Map("rank" -> rankCount, "pvalue" -> 9)) }

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

  val sortedSelectedGenes = additional.toSeq.sortBy(_._2.values).toList.map(_._1)

  // In order to reflect the new found ranking of the genes by the bootstrap, we use the selectedGenes in order of there rank to be written to the selected_nodes file
  PrintFiles.printPattern(config.outputPrefix + ".selected", sortedSelectedGenes, dummyNetworkManager, genePatientMatrix, config.outputFolder, config.inputFolder, additional, malacardsGenes)
  //MutualExclusivityPrintPattern.printPattern(config.outputPrefix + ".selected", rankedGenes.filter { selectedGenes contains _ }, dummyNetworkManager, genePatientMatrix, config.outputFolder, config.inputFolder, additional, malacardsGenes)

  println("Output in " + config.outputFolder + config.outputPrefix + "_networksPPV.bootstrap.tsv")
  println("Network in " + config.outputFolder + config.outputPrefix + ".best_network.html")

}