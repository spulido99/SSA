package be.cmpg.graph

import java.util.LinkedList
import scala.collection.Set

class Node(val gene: Gene) {
  var interactions: Set[(InteractionType.Value, Node)] = _
  var posteriorProbability:Double = _
  var convergenceIteration:Int = -1
  var probabilityHistory:LinkedList[Double] = _
  var score:Double = 1.0
  var bestSubnetwork:(Set[Interaction], Double) = (Set(), 0.0)

  def getInteractions(filter:Option[InteractionType.Value]=None) = {
    if (filter.isEmpty)
      interactions
    else
      interactions.filter(_._1 == filter.get)
  }
  
  override def toString = gene.name + "["+posteriorProbability+"]"
}
