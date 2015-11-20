package be.cmpg.walk.fungus

import be.cmpg.graph.Gene
import be.cmpg.graph.Node
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.walk.SubNetworkSelector

class DFungus(startGene: Gene,
  endGenes: Set[Gene],
  maximumIterations: Int,
  network: NetworkManager[_]) extends SubNetworkSelector(network) {

  val interactionsToVisit = new scala.collection.mutable.HashSet[Interaction]()
  val tree = new SearchTree(startGene, endGenes)
  
  override def getStartGene(): Gene = startGene
  override def getStartNode(): Node = network.getNetwork().getNode(startGene)
  
  def getPossibleInteractions(): List[Interaction] = interactionsToVisit.toList
  
  def getVisitedGenes(): Set[Gene] = tree.getVisitedGenes.toSet
  
  override def selectSubNetwork(): Option[Set[Interaction]] = {

    interactionsToVisit ++= network.getOutgoingInteractionsFor(startGene)

    var i = 0
    while (!interactionsToVisit.isEmpty && i < maximumIterations) {
      val nextInteraction = interactionsToVisit.head
      interactionsToVisit.remove(nextInteraction)

      if (tree.canTakeInteraction(nextInteraction)) {
        nextInteraction.genes
          .filter(!tree.isGeneVisited(_))
          .foreach(gene => interactionsToVisit ++= network.getOutgoingInteractionsFor(gene))

        tree.expand(nextInteraction)

      }
      i += 1
    }

    val result = tree.trim()

    if (result.isEmpty) None else Some(result)
  }
}