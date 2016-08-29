package be.cmpg.cancer.ME

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
import be.cmpg.utils.PatternsManager
import be.cmpg.cancer.Config
import scala.io.Source
import be.cmpg.cancer.ArgumentsParser
import be.cmpg.cancer.DataLoader
import be.cmpg.cancer.PrintFiles
import be.cmpg.walk.WalkerFactory

object MutualExclusivityAnalysis extends App {
  
  val parser = ArgumentsParser.getBasicArgParser("SSA.ME")
  
  // parser.parse returns Option[C]

  parser.parse(args, Config(

    /*
     * Default Values
     */
    iterations = 5000,
    reinforcement = 0.005,
    forgetfulness = 0.995
    
    )) match {

    case Some(config) =>

      val (interactions, network, genePatientMatrix, all_samples, geneList) = DataLoader.loadData(config,false)
      
      val mutualExclusivePatternsManager = new PatternsManager(network)
      
      println("Samples: " + all_samples.size)
      
      val networkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical
        )

      
      val walkers = WalkerFactory.buildHi2iiFungusWalker(geneList.view.toSet, networkManager)
      

      if (config.debug.isDefined)
        networkManager.debug = Some(config.debug.get.map { Gene(_) }.toSet, 20)
      println("Subnetwork selectors: " + walkers.size)

      networkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers),mutualExclusivityPatternManager=Option(mutualExclusivePatternsManager),patterns=config.patterns)

      val rankedGenes = if (config.outputGenes == 0)
        networkManager.rankedGenes.toList
      else
        networkManager.getRankedAllGenes().take(config.outputGenes).toList

      val genesToReport = rankedGenes

      println("\nGenes printed in network: " + rankedGenes.size + " " + rankedGenes.map(_.name))
      
      val malacardsGenes = if (config.otherGeneList==""){Set[Gene]()} else {Source.fromFile(config.otherGeneList).getLines.map(line => Gene(line.split("\t")(0))).toSet}
      
      // Print the resulting files also taking into account malacard genes
      PrintFiles.printPattern(config.outputPrefix, rankedGenes, networkManager, genePatientMatrix, config.outputFolder, config.inputFolder, otherPositiveGeneSetLists=malacardsGenes)

      if (config.patterns){
      println("\n calculating 5 best subnetworks per gene")
      
      mutualExclusivePatternsManager.returnNbestSubnetworks(rankedGenes.toSet,"5bestPatterns",config.outputFolder)
      }
      
      /**
       * P-value calculation
       */
      //val pvalueCalculator = new PValueCalculator(config.pvalIterations, genePatientMatrix, geneList)

    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }

}