package be.cmpg.walk

import be.cmpg.graph.Interaction
import be.cmpg.graph.Gene


class Path(val startGene: Gene) {

  var visitedInteractions = new scala.collection.mutable.MutableList[Interaction]()
  var visitedGenes = new scala.collection.mutable.HashSet[Gene]()
  visitedGenes += startGene

  var currentEndpoint: Gene = startGene

  def getVisitedInteractions() = visitedInteractions.view

  def canTakeInteraction(interaction: Interaction) = {
    interaction.genes.size > 1 && interaction.genes.contains(currentEndpoint) && visitedGenes.intersect(interaction.genes).size == 1
  }
  
  def size() = visitedInteractions.size

  def expand(interaction: Interaction) = {
    if (canTakeInteraction(interaction)) {
      visitedInteractions += interaction
      if (visitedGenes.contains(interaction.from)) {
        visitedGenes.add(interaction.to)
        currentEndpoint = interaction.to
      } else {
        visitedGenes.add(interaction.from)
        currentEndpoint = interaction.from
      }
    } else throw new IllegalArgumentException("Cannot expand randomwalk with the provided interaction.")
  }

  def restart = {
    currentEndpoint = startGene
  }
  
  def reset = {
    currentEndpoint = startGene
    visitedGenes = new scala.collection.mutable.HashSet[Gene]()
    visitedGenes += startGene
    visitedInteractions = new scala.collection.mutable.MutableList[Interaction]()
  }
}