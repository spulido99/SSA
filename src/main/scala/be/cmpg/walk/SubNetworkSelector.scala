package be.cmpg.walk

import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.graph.Node
import scala.util.Random
import be.cmpg.graph.Gene
import scala.collection.Set

/**
 * Represents a method which builds a set of Interaction's given the InteractionNetwork
 *
 */
abstract class SubNetworkSelector(networkManager: NetworkManager[_]) {
  val random = new Random(System.nanoTime())

  /**
   * Select a Set[Interaction] from the supplied NetworkManager.
   * 
   * !!! Maybe we should refactor this to always send a Set but an empty one if the constraints are not met. 
   * 
   * @returns Some(Set[Interaction]) if a valid result is found, returns None() if no valid subnetwork is found.
   */
  def selectSubNetwork(): Option[Set[Interaction]]
  
  def getPossibleInteractions(): List[Interaction]
  
  def getStartGene(): Gene

  def getStartNode(): Node
  
  def getVisitedGenes(): Set[Gene]
  
  def setGeneNumberVariable(value:Double) = {}
  
  def getGeneNumberVariable() = 0.0
}