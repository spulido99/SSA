package be.cmpg.cancer.functional

import be.cmpg.graph.Gene
import scala.collection.JavaConversions._
import be.cmpg.utils.PatternsManager
import be.cmpg.cancer.Config
import be.cmpg.cancer.ME.MutualExclusivityNetworkManager
import be.cmpg.cancer.ArgumentsParser
import be.cmpg.cancer.DataLoader
import be.cmpg.cancer.PrintFiles
import be.cmpg.walk.WalkerFactory

object FunctionalMutualExclusivityAnalysis extends App {
  
  val parser = ArgumentsParser.getBasicArgParser("SSA.FME")

  // parser.parse returns Option[C]

  parser.parse(args, Config(

    /*
     * Default Values
     */
    iterations = 50000,
    reinforcement = 0.0005,
    forgetfulness = 0.9996

    )) match {

    case Some(config) =>

      val (interactions, network, genePatientMatrix, all_samples, geneList) = DataLoader.loadData(config,true)
      
      val mutualExclusivePatternsManager = new PatternsManager(network)
      
      println("Samples: " + all_samples.size)
      
      val networkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical)

      
      val walkers = WalkerFactory.buildHi2iiFungusWalker(geneList.view.toSet, networkManager)
      

      if (config.debug.isDefined)
        networkManager.debug = Some(config.debug.get.map { Gene(_) }.toSet, 20)
      println("Subnetwork selectors: " + walkers.size)

      networkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers), mutualExclusivityPatternManager=Option(mutualExclusivePatternsManager))

      //val rankedGenes = networkManager.getRankedAllGenes().filter { g => networkManager.getPosteriorProbability(g) > 0.95 } .toList
      val rankedGenes = if (config.outputGenes == 0)
        networkManager.rankedGenes.toList
      else
        networkManager.getRankedAllGenes().take(config.outputGenes).toList


      println("\nGenes printed in network: " + rankedGenes.size + " " + rankedGenes.map(_.name))

      PrintFiles.printPattern(config.outputPrefix, rankedGenes, networkManager, genePatientMatrix, config.outputFolder, config.inputFolder)

      /**
       * P-value calculation
       * TO BE IMPLEMENTED
       */
      //val pvalueCalculator = new PValueCalculator(config.pvalIterations, genePatientMatrix, geneList)

    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }

}