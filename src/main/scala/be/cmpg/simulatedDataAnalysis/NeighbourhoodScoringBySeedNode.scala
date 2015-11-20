package be.cmpg.simulatedDataAnalysis

import java.nio.file.Paths
import be.cmpg.graph.Gene
import be.cmpg.graph.NetworkReader
import be.cmpg.graph.Network
import collection.immutable.ListMap
import java.io.FileWriter
import be.cmpg.walk.neighbourhoodScoring.NeighbourhoodTreeGenerator
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.expression.ExpressionNetworkManager

class NeighbourhoodScoringBySeedNode(networkManager:ExpressionNetworkManager, depth:Int = 2) {

  var rootNodeScoreMap = new scala.collection.mutable.HashMap[String, Double]()

  networkManager.getGenes.foreach(mutatedGene => {
    val subNetworkScoreAndGenesInSubNetwork = new NeighbourhoodTreeGenerator(mutatedGene, depth, networkManager).expand
    rootNodeScoreMap.put(mutatedGene.name, subNetworkScoreAndGenesInSubNetwork._1)
  })

  val orderedRootNodeScoreMap = ListMap(rootNodeScoreMap.toList.sortBy { _._2 }: _*)

  val selectedRootNodesWithValue = orderedRootNodeScoreMap.map(input => input._1 + "\t" + input._2.toString).toList

  val rankedGenes = selectedRootNodesWithValue.map(element => element.split("\t")(0)).toList
  
}