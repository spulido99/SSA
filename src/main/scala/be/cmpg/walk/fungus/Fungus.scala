package be.cmpg.walk.fungus

import scala.collection.mutable.HashSet
import scala.util.Random
import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.walk.Path
import be.cmpg.graph.Network
import be.cmpg.utils.StatUtils
import be.cmpg.graph.Interaction
import scala.collection.Set
import be.cmpg.graph.Node
import be.cmpg.graph.InteractionType
import collection.JavaConversions._
import be.cmpg.graph.Interaction

class Fungus(startGene: Gene,
    endGenes: Set[Gene] = Set(),
    _geneNumberVariable: Double = 3,
    network: NetworkManager[_]) extends SubNetworkSelector(network) {

  if (!network.getGenes.contains(startGene))
    throw new IllegalArgumentException("Starting gene " + startGene + " does not exist in the network");

  val startNode = network.getNetwork().getNode(startGene)
  var visitedNodes = new HashSet[Node]
  private var extendableNodes = new HashSet[Node]
  //val posibleInteractions = new HashSet[(Node, Node)]

  var currentScore: Double = _
  var geneNumberVariable = _geneNumberVariable

  override def setGeneNumberVariable(value: Double) = {
    geneNumberVariable = value
  }

  override def getGeneNumberVariable() = geneNumberVariable

  override def getStartGene(): Gene = startGene

  override def getStartNode(): Node = startNode

  override def getPossibleInteractions(): List[Interaction] = throw new RuntimeException

  override def getVisitedGenes() = visitedNodes.map(_.gene)

  def allowAddMoreNodes(): Boolean = {
    currentScore < geneNumberVariable && !extendableNodes.isEmpty
  }

  // Something is wrong here... The interactions are too random (they do not form subnetworks). Might have something to do with
  def getRandomInteraction(): Option[(Node, Node)] = {
    // Note that every interaction is viewed as an undirected interaction
    val randomNode = extendableNodes.toList(Random.nextInt(extendableNodes.size))
    val posibleInteractions = randomNode.getInteractions().filter(e => !visitedNodes.contains(e._2))

    if (posibleInteractions.isEmpty) {

      extendableNodes -= randomNode
    }

    //   val sum = posibleInteractions.foldLeft(0.0)((x, y) => x + math.max(y._2.posteriorProbability, 0.01))
    val sum = posibleInteractions.foldLeft(0.0)((x, y) => x + y._2.posteriorProbability)
    var ran = Random.nextDouble() * sum
    var selected: (InteractionType.Value, Node) = null

    for (i <- posibleInteractions) {
      ran -= i._2.posteriorProbability
      //      ran -= math.max(i._2.posteriorProbability, 0.01)
      selected = i
      if (ran < 0) {
        return Some(randomNode, selected._2)
      }
    }
    None
  }

  def selectSubNetwork(): Option[Set[Interaction]] = {

    /*
     * Initialize used variables
     */
    currentScore = startNode.score
    val visitedInteractions = new HashSet[(Node, Node)]

    visitedNodes = new HashSet[Node]
    visitedNodes += startNode

    extendableNodes = new HashSet[Node]
    extendableNodes += startNode

    //posibleInteractions.clear()
    //posibleInteractions ++= startNode.getInteractions().map(e => (startNode, e._2))

    /*
     * Select subnetwork
     */

    // Done and Nonecounter are used to terminate the process when no additional interaction can be added to the network after multiple tries
    var done = false
    var NoneCounter = 0

    while (allowAddMoreNodes && !(done)) {
      val interaction = getRandomInteraction

      if (interaction.isDefined) {
        NoneCounter = 0
        val selectedNode = if (visitedNodes contains interaction.get._1)
          interaction.get._2
        else
          interaction.get._1

        //posibleInteractions -= interaction.get

        visitedNodes += selectedNode
        extendableNodes += selectedNode
        visitedInteractions += interaction.get

        //posibleInteractions ++= selectedNode.getInteractions()
        //  .map(e => (selectedNode, e._2))
        //  .filter(e => !visitedNodes.contains(interaction.get._1) || !visitedNodes.contains(interaction.get._2))

        currentScore += selectedNode.score
      } else {
        NoneCounter = NoneCounter + 1
        if (NoneCounter > 50) { done = true }
      }
    }

    prune()
    Some(visitedInteractions.map(i => Interaction(i._1.gene, i._2.gene)).toSet)
  }

  private def prune() = {
    /*
     * This function should prune the paths that do not take to an end gene
     */
    if (!endGenes.isEmpty) {
      // TODO do pruning
    }
  }

  private def connectedSubnetwork(interactions: Set[Interaction]): Boolean = {
    var connected = true
    interactions.foreach(interaction => {
      val otherInteractions = interactions.-(interaction)
      val otherGenes = otherInteractions.toList.map(interaction => interaction.genes).flatten.map(gene => gene.name)
      if (!(otherGenes.contains(interaction.from.name)) && !(otherGenes.contains(interaction.to.name))) { connected = false }
    })
    connected
  }

  /*
  def cycle(): Set[Interaction] = {
    val genes = new scala.collection.mutable.HashSet[Gene]()

    network.genes.foreach(genes.add(_))

    val tree = new SearchTree(startGene, endGenes)
    
    var currentGene = startGene

    while(!genes.isEmpty && score < maximumScore){
      
      network.getOutgoingInteractions(currentGene)
    }
    null
  }
  * 
  */
}