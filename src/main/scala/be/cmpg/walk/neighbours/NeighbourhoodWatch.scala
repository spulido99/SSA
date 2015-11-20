package be.cmpg.walk.neighbours

import be.cmpg.graph.Network
import be.cmpg.graph.Node
import be.cmpg.graph.Gene
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.NetworkManager

class NeighbourhoodWatch(startGene: Gene,
  depth: Int,
  network: NetworkManager[Interaction]) extends SubNetworkSelector(network) {

  val visitedInteractions = new scala.collection.mutable.HashSet[Interaction]()
  val visitedGenes = new scala.collection.mutable.HashSet[Gene]()
  
  def getPossibleInteractions(): List[Interaction] = List.empty
  
  def getVisitedGenes(): Set[Gene] = visitedGenes.toSet
  
  override def getStartGene(): Gene = startGene
  override def getStartNode(): Node = network.getNetwork().getNode(startGene)
  
  override def selectSubNetwork(): Option[Set[Interaction]] = {

    var nextInteractionsToVisit = new scala.collection.mutable.Stack[Interaction]()

    network.getOutgoingInteractionsFor(startGene).foreach(nextInteractionsToVisit.push(_))
    visitedGenes += startGene

    var i = 0

    while (i < depth) {
      val newInteractionsToVisit = new scala.collection.mutable.Stack[Interaction]()

      while (!nextInteractionsToVisit.isEmpty) {
        val currentInteraction = nextInteractionsToVisit.pop
        val currentEndGene =
          if (visitedGenes.contains(currentInteraction.from) && visitedGenes.contains(currentInteraction.to)) None
          else if (visitedGenes.contains(currentInteraction.from)) Some(currentInteraction.to)
          else if (visitedGenes.contains(currentInteraction.to)) Some(currentInteraction.from)
          else throw new IllegalStateException("We can never access an non-linked interaction.")

        if (currentEndGene.isDefined) {
          val nextInteractions = network.getOutgoingInteractionsFor(currentEndGene.get).filter(interaction => !visitedInteractions.contains(interaction))

          nextInteractions.foreach(x => {
            newInteractionsToVisit.push(x)
            visitedInteractions.add(x)
          })
          visitedGenes += currentEndGene.get
        }
      }

      nextInteractionsToVisit = newInteractionsToVisit
      i += 1
    }

    if (visitedInteractions.isEmpty) None
    else Some(visitedInteractions.toSet)
  }

}