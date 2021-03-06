package be.cmpg.cancer

import be.cmpg.graph.Network
import scala.collection.Set
import be.cmpg.graph.Gene
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Interaction
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.walk.fungus.Fungus

class PancancerNetworkManager(network: Network,
                              genePValues: Map[Gene, PValueInfo],
                              pheromone: Double,
                              evaporation: Double,
                              ranked: Boolean,
                              convergenceThreshold:Double) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, ranked: Boolean, convergenceThreshold:Double) {

  def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    val nodes = subnetwork.map(_.genes).flatten
    val startGene = if (selector.isDefined) selector.get.getStartGene() else Gene("")
    
    if (selector.isDefined && nodes.size < selector.get.getGeneNumberVariable) return Double.NaN
    
    var result = 0.0
    var resultMinusStart = 0.0
    nodes.foreach{ node => {
         val pValues = genePValues.get(node)
         val value = if (pValues.isDefined) pValues.get.max else 0.0
         result += value
         if (node != startGene) {
        	 resultMinusStart += value
         }
      }}
      
    result /= nodes.size
    resultMinusStart /= nodes.size - 1
    // Return NaN if the network would be better off without the starting gene
    if (resultMinusStart < result) Double.NaN else result
  }

}