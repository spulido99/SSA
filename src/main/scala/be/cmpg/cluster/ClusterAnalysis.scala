package be.cmpg.cluster

import scala.collection.JavaConversions._
import be.cmpg.utils.PatternsManager
import be.cmpg.walk.WalkerFactory

object MushthofaAnalysis extends App {

  val parser = ArgumentsParser.getBasicArgParser("Cluster")

  // parser.parse returns Option[C]

  parser.parse(args, PvalueConfig()) match {

    case Some(config) =>

      val (interactions, network, pValueMap, geneList) = DataLoader.loadData(config)

      println(geneList.size)

      val clusterPatternsManager = new PatternsManager(network)

      println("number of seed genes: " + geneList.size)

      val networkManager = new ClusterNetworkManager(
        network = network,
        genePValueMatrix = pValueMap,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank)

      val walkers = WalkerFactory.buildHi2iiFungusWalker(geneList.view.toSet, networkManager)
      
      println("Number of Subnetwork selectors: " + walkers.size)

      networkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers), mutualExclusivityPatternManager = Option(clusterPatternsManager), patterns = config.patterns)

      //val rankedGenes = networkManager.getRankedAllGenes().filter { g => networkManager.getPosteriorProbability(g) > 0.95 } .toList
      val rankedGenes = if (config.outputGenes == 0)
        networkManager.rankedGenes.toList
      else
        networkManager.getRankedAllGenes().take(config.outputGenes).toList

      val genesToReport = rankedGenes

      println("\nGenes printed in network: " + rankedGenes.size + " " + rankedGenes.map(_.name))

      PValuesPrintPattern.printPattern(config.outputPrefix, rankedGenes, networkManager, pValueMap, config.outputFolder, config.inputFolder)

      if (config.patterns) {
        println("\n calculating 5 best subnetworks per gene")

        clusterPatternsManager.returnNbestSubnetworks(rankedGenes.toSet, "5bestPatterns", config.outputFolder)
      }

    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }

}