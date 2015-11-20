package be.cmpg.graph.interaction

import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.walk.Path
import be.cmpg.graph.Gene
import java.io.PrintStream
import scala.collection.mutable.HashMap
import scala.collection.Iterable
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.graph.Interaction
import scala.collection.Map
import scala.collection.Set

abstract class NetworkManager[T](network: Network) {

  def getNetwork() = network
  
  def scoreWalker(selector: SubNetworkSelector) : Option[(Set[Interaction], Double)] ={
    val sn = selector.selectSubNetwork
    if (sn.isDefined) Some(sn.get, scoreSubnetwork(sn.get, Some(selector))) else None
  }

  def scoreSubnetwork(subnetwork: Set[Interaction], selector:Option[SubNetworkSelector]=None) : Double
  
  def updateScores(subnetworkScores: Traversable[(Set[Interaction], Double)]) : Map[T, Double]

  def evaporate(scores:Map[T, Double] = Map())

  def getPosteriorProbability(obj: T): Double

  def printResults(writer: PrintStream) =
    writer.println(
      network
        .getNodes()
        .map(x => x.gene.name.toString + "\t" + x.posteriorProbability)
        .mkString("\n"))

  def getGenes = network.genes

  def getOutgoingInteractionsFor(gene: Gene): Set[Interaction] = network.getOutgoingInteractions(gene)
  
  //def setGeneScore(gene: Gene, score: Double) = geneScoreMap.put(gene, score)
  
  def getAllInteractions = network.interactions
  
  //def getGeneScore(gene: Gene) : Double = if (geneScoreMap.contains(gene)) geneScoreMap(gene) else 1.0
  
  def getRandomInteraction(selector: SubNetworkSelector) : Option[Interaction]
  
}