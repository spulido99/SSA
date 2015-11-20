package be.cmpg.walk.fungus

import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction

class SearchTree(start: Gene, endGenes: Set[Gene]) {

  private val visitedInteractions = new scala.collection.mutable.HashSet[Interaction]()
  private val visitedGenes = new scala.collection.mutable.HashSet[Gene]()
  private val leafGenes = new scala.collection.mutable.HashSet[Gene]()

  visitedGenes.add(start)

  def expand(interaction: Interaction) = {
    if (canTakeInteraction(interaction)) {
      visitedInteractions += interaction
      if (visitedGenes.contains(interaction.from)) {
        visitedGenes.add(interaction.to)
        leafGenes.remove(interaction.from)
        leafGenes.add(interaction.to)
      } else {
        visitedGenes.add(interaction.from)
        leafGenes.remove(interaction.to)
        leafGenes.add(interaction.from)
      }
    } else throw new IllegalArgumentException("Could not connect the given interaction to the current search tree.")
  }

  def getLeaves = leafGenes.view

  def getVisitedGenes() = visitedGenes
  
  def isLeafGene(gene: Gene) = leafGenes.contains(gene)

  def canTakeInteraction(interaction: Interaction) = {
    visitedGenes.intersect(interaction.genes).size == 1
  }

  def isGeneVisited(gene: Gene) = visitedGenes.contains(gene)

  def contains(interaction: Interaction) = visitedInteractions.contains(interaction)

  def trim(): Set[Interaction] = {
    val nonEndLeafGenes = leafGenes -- endGenes

    while (!nonEndLeafGenes.isEmpty) {
      val currentLeafGene = nonEndLeafGenes.head

      val interactionToCurrentLeafGene = 
        	visitedInteractions
        		.find(interaction => interaction.genes.contains(currentLeafGene))
        		.get
        		
      val otherEndGene = (interactionToCurrentLeafGene.genes - currentLeafGene).head

      nonEndLeafGenes.remove(currentLeafGene)
      visitedInteractions.remove(interactionToCurrentLeafGene)

      if (otherEndGene != start &&
        !endGenes.contains(otherEndGene) &&
        visitedInteractions.filter(interaction => interaction.genes.contains(otherEndGene)).size == 1) {
        nonEndLeafGenes.add(otherEndGene)
      }
    }

    visitedInteractions.toSet
  }

}