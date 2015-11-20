package be.cmpg.graph

import be.cmpg.graph.interaction.NetworkManager
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import scala.collection.Set

class Network(val interactions: Set[Interaction], geneList: Set[Gene] = null) {

  val size = interactions.size

  lazy val genes = interactions.map(interaction => Set(interaction.from, interaction.to)).flatten[Gene]

  val nodeMap = {

    val nodeMap = genes.map(gene => (gene, new Node(gene))).toMap

    val interactionMapMut = genes.map(gene => (gene, new HashSet[(InteractionType.Value, Node)])).toMap[Gene, HashSet[(InteractionType.Value, Node)]]
    interactions.foreach(interaction => {
      interactionMapMut(interaction.from) += ((InteractionType.Outgoing, nodeMap(interaction.to)))
      interactionMapMut(interaction.to) += ((InteractionType.Incoming, nodeMap(interaction.from)))

      if (interaction.direction == "undirected") {
        interactionMapMut(interaction.from) += ((InteractionType.Incoming, nodeMap(interaction.to)))
        interactionMapMut(interaction.to) += ((InteractionType.Outgoing, nodeMap(interaction.from)))
      }
    })

    genes.foreach(gene => nodeMap(gene).interactions = interactionMapMut(gene).toSet)

    nodeMap
  }

  /**
   * @Deprecated use node.interactions
   */
  @Deprecated
  def getOutgoingInteractions(gene: Gene): Set[Interaction] = {
    val node = nodeMap(gene)
    node.getInteractions(Some(InteractionType.Outgoing)).map( e => Interaction(node.gene, e._2.gene))
  }

  def getNodes() = nodeMap.values
  
  def getNode(gene:Gene) = nodeMap(gene)
  
  @Deprecated
  def getAllInteractions(gene: Gene): Set[Interaction] = {
    val node = nodeMap(gene)
    node.getInteractions().map( e => Interaction(node.gene, e._2.gene))
  }

  @Deprecated
  def getIncommingInteractions(gene: Gene): Set[Interaction] = {
    val node = nodeMap(gene)
    node.getInteractions(Some(InteractionType.Incoming)).map( e => Interaction(node.gene, e._2.gene))
  }

}