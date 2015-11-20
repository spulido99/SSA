package be.cmpg.walk

import be.cmpg.graph.Gene
import be.cmpg.graph.Node
import be.cmpg.graph.Network
import scala.util.Random
import java.util.concurrent.Callable
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.NetworkManager

class Ant[T](restartFrequency: Double,
  maximumIterations: Int,
  startGene: Gene,
  endGenes: Set[Gene],
  network: NetworkManager[_]) extends SubNetworkSelector(network) {

  val path = new Path(startGene)
  
  override def getVisitedGenes() = path.visitedGenes.toSet
  
  override def getStartGene(): Gene = startGene
  
  override def getStartNode(): Node = network.getNetwork().getNode(startGene)
  
  override def selectSubNetwork(): Option[Set[Interaction]] = {

    var iteration = 0

    if (network.getOutgoingInteractionsFor(startGene).filter(interaction => path.canTakeInteraction(interaction)).isEmpty) return None

    def checkIfPathIsOk(path: Path) = {
      if (path.size == 0) false
      else endGenes.contains(path.currentEndpoint) && path.getVisitedInteractions.last.regulatory == "regulatory"
    }

    while (!checkIfPathIsOk(path) && iteration < maximumIterations) {
      if (Random.nextDouble <= restartFrequency) {
        path.reset
      } else {
        val next = network.getRandomInteraction(this)
        if (next.isDefined) path.expand(next.get) else path.reset
      }
      iteration += 1
    }

    if (checkIfPathIsOk(path)) Some(path.getVisitedInteractions.toSet) else None
  }

  override def getPossibleInteractions(): List[Interaction] = {
    network.getOutgoingInteractionsFor(path.currentEndpoint).filter(interaction => path.canTakeInteraction(interaction)).toList
  }

}