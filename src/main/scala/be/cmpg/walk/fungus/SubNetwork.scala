package be.cmpg.walk.fungus

import scala.collection.mutable.HashSet
import be.cmpg.graph.Interaction
import be.cmpg.graph.Interaction
import be.cmpg.graph.Gene

class SubNetwork(val startGene: Gene) {

  var visitedInteractions = new scala.collection.mutable.MutableList[Interaction]()
  var visitedGenes = new scala.collection.mutable.HashSet[Gene]()
  visitedGenes += startGene
  
  def canTakeInteraction(interaction: Interaction) = {
    visitedGenes.intersect(interaction.genes).size > 0
  }
  
  def getVisitedInteractions() = visitedInteractions
  
  def size() = visitedInteractions.size

  def expand(interaction: Interaction) = {
    if (canTakeInteraction(interaction)) {
      
      visitedInteractions += interaction
      
      if (visitedGenes.contains(interaction.from)) {
        visitedGenes.add(interaction.to)
      } else {
        visitedGenes.add(interaction.from)
      }
    } else throw new IllegalArgumentException("Cannot expand randomwalk with the provided interaction.")
  }

}