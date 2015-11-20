package be.cmpg.cancer

import scopt.OptionParser
import java.io.File
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

object MutualExclusivityAnalysis extends App {

  case class Config(
    iterations: Int,
    reinforcement: Double,
    forgetfulness: Double,
    refNetwork: Seq[String],
    useRank: Boolean,
    hypergeometricTest: Boolean,
    minMutPerGene: Int,
    seedGenesMutations: Int,
    outputPrefix: String,
    input: Option[String] = None,
    outputGenes:Option[Int] = None,
    subtype:Option[String] = None,
    subtypeFile:Option[File] = None,
    debug:Option[Seq[String]] = None)

  val parser = new scopt.OptionParser[Config]("SSA.ME.") {
    head("Mutual exclusivity analysis by Small Subnetwork Analysis", "0.1")
    opt[Int]('i', "iterations") action { (x, c) =>
      c.copy(iterations = x)
    } text ("iterations is the number of iterations (default: 5000)")

    opt[Double]('r', "reinforcement") action { (x, c) =>
      c.copy(reinforcement = x)
    } text ("The reinforncement for succesful learning (default: 0.005)")

    opt[Double]('f', "forgetfulness") action { (x, c) =>
      c.copy(forgetfulness = x)
    } text ("The forgetfulness after each iteration (default: 0.994)")

    opt[Int]('s', "seedGenesMutations") action { (x, c) =>
      c.copy(seedGenesMutations = x)
    } text ("The number of mutated samples required in a gene to be a seed (default: 1)")

    opt[String]('o', "outputPrefix") action { (x, c) =>
      c.copy(outputPrefix = x)
    } text ("The prefix for the analisis output files (default: ME)")

    opt[Seq[String]]('n', "refNetwork") action { (x, c) =>
      c.copy(refNetwork = x)
    } text ("The reference network (HT, hiII14, ppi, cell-map, nci-nature, reactome or/and files. HT, hiII14 and reactoe used by default.)")

    /*    
    opt[Int]('g', "maxOutputGenes") action { (x, c) =>
      c.copy(maxOutputGenes = x)
    } text ("Jump by ")
    */

    opt[Boolean]('u', "useRank") action { (x, c) =>
      c.copy(useRank = x)
    } text ("Rank the gene scores instead of the scaled values (default: true)")
    
    opt[Boolean]('h', "hypergeometricTest") action { (x, c) =>
    c.copy(hypergeometricTest = x)
    } text ("Use an hypergeometric test to search for mutual exclusivity. Will not find low mutations genes. (default: false)")
    
    opt[Int]('g', "minMutPerGene") action { (x, c) =>
      c.copy(minMutPerGene = x)
    } text ("Ignore genes mutated in less than 'g' samples (default: 1)")

    opt[String]('m', "input_matrix") required () action { (x, c) =>
      var inputMatrix = x
      if (!x.endsWith(".m2"))
        inputMatrix = x + ".m2"
      c.copy(input = Some(inputMatrix))
    } text ("Mutation matrix .m2 file (created with option ME_input )")

    opt[Int]('p', "outputGenes") action { (x, c) =>
      c.copy(outputGenes = Some(x))
    } text ("Number of nodes to be published to the html network output (default: 50)")
    
    opt[File]("subtypeFile") action { (x, c) =>
      c.copy(subtypeFile = Some(x))
    } text ("Samples per subtype (Tab Delimited File with 2 columns: Sample Subtype")
    
    opt[String]("subtype") action { (x, c) =>
      c.copy(subtype = Some(x))
    } text ("subtype to analyse")
    
    opt[Seq[String]]('d', "debug") action { (x, c) =>
      c.copy(debug = Some(x))
    } text ("Debug mode. List of genes that want to be observed (e.g. -d TP53,MYC,")
    
    help("help")

  }
  // parser.parse returns Option[C]

  val helper = new CancerHelper
  
  parser.parse(args, Config(
      
      /*
       * Default Values
       */
      iterations=50000,
      reinforcement=0.0005,
      forgetfulness=0.9996,
      refNetwork=Seq("HT", "hiII14", "reactome"), 
      useRank=true, 
      hypergeometricTest=false, 
      minMutPerGene=3, 
      seedGenesMutations=1, 
      outputPrefix="ME")) match {
    
    case Some(config) =>
      
      val interactions = helper.loadNetwork(config.refNetwork)
      val network = new Network(interactions)
      val translateGenesToEntrez: Map[String, String] = network.genes.map(g => (g.name, g.name)).toMap

      println("Interactions: "+interactions.size)
      
      val genePatientMatrix = {

        println("Loading mutation matrix file...")
        val genePatientMatrix = helper.loadMutationMatrixFiles(config.input.get)
        
        if (config.subtype.isDefined && config.subtypeFile.isDefined) {
          
          //subtypeFile = TCGA_BRCA_clin_PAM50.txt
          val samplesInfo = new CSVReader(new FileReader(config.subtypeFile.get), '\t')
                        .readAll()
                        .map {fields => SampleInfo(fields(0).trim.replaceAll("\\.", "-"), fields(fields.length - 1).trim) }
                        .filter { _.subtype.equalsIgnoreCase(config.subtype.get) }
                        .map { _.sample }
          
                                                
          genePatientMatrix.retain{ (pk, p) => samplesInfo.contains(pk.sample) } 
        }
        
        genePatientMatrix.toMap
      }

      val all_samples = genePatientMatrix.map(_._1.sample).toSet
      println("Samples: "+all_samples.size)
      /*
      println("===========================")
      all_samples.foreach { sample => println(sample) }
      println("===========================")
      */
      
      val geneList = new HashSet[Gene]
      network.genes.foreach(gene => {
        var count = 0
        for (sample <- all_samples) {
          if (genePatientMatrix contains (PolimorphismKey(gene, sample))) {
            count += 1
          }
        }

        if (count >= config.seedGenesMutations) {
          geneList.add(gene)
        }
      })

      val networkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        hypergeometricTest = config.hypergeometricTest)

      val walkers: Set[SubNetworkSelector] = geneList.map(gene =>
        //new Fungus(
        new Hi2iiFungus(
          startGene = gene,
          geneNumberVariable = 2 + StatUtils.getRandomPoisson(1),
          //geneNumberVariable = 3 + StatUtils.getRandomPoisson(0.5), // 0 - 60%, 1 - 30%, 2 - 10%
          //geneNumberVariable = 3,
          network = networkManager)).toSet

      if (config.debug.isDefined)
        networkManager.debug = Some(config.debug.get.map { Gene(_) }.toSet, 20)
      //networkManager.debug = Some(Set(Gene("TP53")), 20)
      println("Subnetwork selectors: "+walkers.size)
          
      networkManager.run(config.iterations, geneList.view.toSet, Runtime.getRuntime().availableProcessors() / 2, defwalkers = Some(walkers))

      //val rankedGenes = networkManager.getRankedAllGenes().filter { g => networkManager.getPosteriorProbability(g) > 0.95 } .toList
      val rankedGenes = if (config.outputGenes.isEmpty) 
                          networkManager.rankedGenes.toList 
                        else 
                          networkManager.getRankedAllGenes().take(config.outputGenes.get).toList

      val genesToReport = rankedGenes
                          /*.map { g => network.getAllInteractions(g).map { i => i.genes }.flatten }.flatten
                          .filter { g => {
                            rankedGenes.contains(g) || network.getAllInteractions(g).map { i => i.genes }.flatten.toSet.intersect(rankedGenes).size > 2
                          }}*/
      
      /*val genesToReport = rankedGenes.toList
                          .map { g => network.getAllInteractions(g).map { i => i.genes }.flatten }.flatten
                          .groupBy(identity).mapValues(_.size)
                          .filter(e => rankedGenes.contains(e._1) || e._2 > 2)
                          .map(_._1)
                          .toSet*/
      
      
      println("Genes printed in network: "+rankedGenes.size + " "+rankedGenes.map(_.name))
      
      MutualExclusivityPrintPattern.printPattern(config.outputPrefix, rankedGenes, networkManager, genePatientMatrix)
      
    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }

}