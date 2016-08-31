package be.cmpg.cancer.pancancer

import be.cmpg.graph.Network
import scala.collection.Set
import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Interaction
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.utils.weightByFlatInitialProbability

class PancancerNetworkManager(network: Network,
                              genePValues: Map[Gene, PValueInfo],
                              pheromone: Double,
                              evaporation: Double,
                              ranked: Boolean,
                              initialProb:Double,
                              convergenceThreshold:Double) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, new weightByFlatInitialProbability(network,initialProb), ranked: Boolean, convergenceThreshold:Double) {

  def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    val genes = subnetwork.map(_.genes).flatten
    val startGene = if (selector.isDefined) selector.get.getStartGene() else Gene("")
    
    if (selector.isDefined && genes.size < selector.get.getGeneNumberVariable) return Double.NaN
    
    var result = 0.0
    var resultMinusStart = 0.0
    genes.foreach{ gene => {
         val pValues = genePValues.get(gene)
         val value = if (pValues.isDefined) pValues.get.max else 0.0
         result += value
         if (gene != startGene) {
        	 resultMinusStart += value
         }
      }}
      
    result /= genes.size
    resultMinusStart /= genes.size - 1
    // Return NaN if the network would be better off without the starting gene
    if (resultMinusStart < result) Double.NaN else result
  }

}