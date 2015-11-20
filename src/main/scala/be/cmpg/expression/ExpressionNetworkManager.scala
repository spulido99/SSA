package be.cmpg.expression

import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import scala.collection.Set
import be.cmpg.graph.Gene
import be.cmpg.walk.SubNetworkSelector

class ExpressionNetworkManager(network: Network,
  pheromone: Double = 0.05,
  evaporation: Double = 0.996) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double) {

  var geneExpression: Map[Gene, Double] = _
  
  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    var score = 0.0;
    subnetwork.map(_.genes).flatten[Gene].foreach(gene => score += geneExpression(gene).abs)
    score
  }

}