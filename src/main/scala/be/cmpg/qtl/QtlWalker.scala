package be.cmpg.qtl

import be.cmpg.graph.Gene
import be.cmpg.graph.Node
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.walk.Path
import be.cmpg.graph.Interaction
import scala.collection.mutable.LinkedList
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ListBuffer
import scala.collection.mutable.TreeSet
import scala.util.Random
import scala.collection.mutable.HashSet
import be.cmpg.graph.Node

/*
 * The selected subnetwork is None if the start or end Gene does not have any marker
 * that affect it
 * 
 * 
 * FIXME What do we do if there are no other markers around the gene marker?
 */
class QtlWalker(startGene: Gene, 
  qtlHalfRegionSize: Int, 
  pathSize: Int, 
  network: NetworkManager[_]) extends SubNetworkSelector(network) {

  
  val path = new Path(startGene)
  val startNode = network.getNetwork().getNode(startGene)
  val interactions = new HashSet[Interaction]

  override def getStartGene(): Gene = startGene
  override def getStartNode(): Node = startNode
  
  override def getVisitedGenes() = interactions.map(_.genes).flatten

  override def getPossibleInteractions(): List[Interaction] = {
    network.getOutgoingInteractionsFor(path.currentEndpoint).filter(interaction => path.canTakeInteraction(interaction)).toList
  }

  def selectSubNetwork(): Option[Set[Interaction]] = {

    val geneOfInterest = startGene
    val (startGene2Markers, startMarkers) = getFlankingMarkers(geneOfInterest)

    if (startGene2Markers.isEmpty)
      None
    
    startGene2Markers.foreach(interactions += _)
    startMarkers.sliding(2).foreach(mt => interactions += Interaction(mt(0), mt(1), QtlConstants.startLoci))

    var steps = 0
    while (steps < pathSize) {
      
      val next = network.getRandomInteraction(this)
      if (next.isDefined)
        path.expand(next.get)
      steps += 1
    }
    
    if (startMarkers.contains(path.currentEndpoint))
      None // return None if the random path returned to the original sets of markers

    val (endGene2Markers, endMarkers) = getFlankingMarkers(path.currentEndpoint)

    if (endGene2Markers.isEmpty)
      None
    
    endGene2Markers.foreach(interactions += _)
    endMarkers.sliding(2).foreach(mt => interactions += Interaction(mt(0), mt(1), QtlConstants.endLoci))
    
    Some(interactions.view.toSet)
  }
  
  private def getFlankingMarkers(geneOfInterest: be.cmpg.graph.Gene): (List[be.cmpg.graph.Interaction], List[be.cmpg.graph.Gene]) = {
    val gene2markers = network.getOutgoingInteractionsFor(geneOfInterest).filter(_.typ == "marker").toList.sortBy(_.to.name)

    val chromosome = gene2markers(0).to.get(QtlConstants.chromosome).get
    val startPosition = chromosome + math.max(0, gene2markers(0).to.get(QtlConstants.position).get.toInt - qtlHalfRegionSize)
    val endPosition = chromosome + (gene2markers(0).to.get(QtlConstants.position).get.toInt + qtlHalfRegionSize)

    val markers = network.getGenes.filter(node => startPosition <= node.name && node.name <= endPosition).toList.sortBy(_.name)
    
    (gene2markers, markers)
  }

}