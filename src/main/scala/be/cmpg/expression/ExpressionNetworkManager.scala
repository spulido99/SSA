package be.cmpg.expression

import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import scala.collection.Set
import be.cmpg.graph.Gene
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.utils.NetworkWeighting

class ExpressionNetworkManager(network: Network,
  initialProb:Double =0.5,
  pheromone: Double = 0.05,
  evaporation: Double = 0.996,
  weightingScheme: NetworkWeighting) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double,weightingScheme: NetworkWeighting ) {

  var geneExpression: Map[Gene, Double] = _
  
  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    var score = 0.0;
    subnetwork.map(_.genes).flatten[Gene].foreach(gene => score += geneExpression(gene).abs)
    score
  }

}