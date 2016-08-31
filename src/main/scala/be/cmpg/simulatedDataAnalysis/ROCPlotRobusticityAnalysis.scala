package be.cmpg.simulatedDataAnalysis

import be.cmpg.graph.NetworkReader
import be.cmpg.graph.Network
import java.nio.file.Paths
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Gene
import scala.collection.mutable.HashSet
import be.cmpg.graph.Interaction
import be.cmpg.graph.Interaction
import java.io.FileWriter
import be.cmpg.expression.ExpressionNetworkManager
import be.cmpg.utils.weightByFlatInitialProbability

object ROCPlotRobusticityAnalysis extends App {
  val repeats = 100

  val threadpool = java.util.concurrent.Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1)
  //val threadpool = java.util.concurrent.Executors.newSingleThreadExecutor()

  val interactions = NetworkReader.fromFile("src/main/resources/be/cmpg/graph/network_1.txt")
  val allGenes = interactions.map(interaction => Set(interaction.from, interaction.to)).flatten[Gene]
  val size = 10

  def changedEdgesNetworkManager(removeChance: Double, base_networkManager: NodeCostNetworkManager) : ExpressionNetworkManager = {

    val newInteractions = new HashSet[Interaction]

    val geneList = allGenes.toArray

    for (i <- interactions) {
      if (math.random < removeChance) {
        val from = geneList((math.random * geneList.size).toInt)
        val to = geneList((math.random * geneList.size).toInt)
        newInteractions += Interaction(from, to)
      } else {
        newInteractions += i
      }
    }
    
    val network = new Network(newInteractions.view.toSet)
    
    val networkManager = new ExpressionNetworkManager(network = network,weightingScheme = new weightByFlatInitialProbability(network,0.5))
    networkManager.getNetwork.getNodes.foreach(node => node.score = base_networkManager.getNetwork.getNode(node.gene).score)
    
    networkManager
  }

  val base_network = new Network(interactions)
  

  val writer = new FileWriter("Robusticity_" + size + "_output.txt")

  writer.write("Rep\tN\tNoise\taspecificity\tsensitivity")

  val changeLevels = List(0.0, 0.1, 0.25, 0.5)

  for (r <- 0 to repeats) {

    println("---> " + r)

    val base_networkManager = new ExpressionNetworkManager(network = base_network,weightingScheme = new weightByFlatInitialProbability(base_network,0.5))
    val (importantGenesR, geneExpression) = ROCPlotSimulatedData.simulateGeneScores(base_networkManager, size)

    val importantGenes = importantGenesR.map(_.name).toSet
    base_networkManager.geneExpression = geneExpression.view.toMap;
    
    println("- Network Manager created")

    for (changeLevel <- changeLevels) {

      val networkManager = changedEdgesNetworkManager(changeLevel, base_networkManager)
      networkManager.geneExpression = base_networkManager.geneExpression;
      val allGenes = networkManager.getGenes.map(_.name).toSet

      print("- Start rank by Phac. Perturbation Level: "+changeLevel)
      val rankedGeneByPhac = ROCPlotSimulatedData.getRankedGenesByPhac(networkManager, threadpool = threadpool).toList.map(_.name);
      println("Done!")

      val phacSelected = new HashSet[String]
      for (n <- 0 to rankedGeneByPhac.size - 1) {

        phacSelected += rankedGeneByPhac(n)
        val phacROCValues = new ROCPlotValuesCreator(phacSelected.view.toSet, importantGenes, allGenes).get

        writer.write("\n" + r + "\t" + n + "\t" + changeLevel + "\t" + phacROCValues.sensitivity + "\t" + phacROCValues.aspecificity)
      }
    }
  }

}