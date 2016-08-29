package be.cmpg.cluster

import be.cmpg.graph.Network
import scala.io.Source
import be.cmpg.graph.Gene
import java.io.File
import scala.collection.mutable.HashSet

object DataLoader {
  
    def loadData(config: PvalueConfig) = {
    val interactions = NetworkFactory.loadNetwork(new File(config.refNetwork))
    // Set up the network, taking into account the excluded genes
    val excludedGenes = if (!(config.excludedGenes =="")) {Source.fromFile(new File(config.excludedGenes)).getLines} else {List()}
    val network = new Network(interactions.filter(interaction => !(config.excludedGenes.contains(interaction.from.name) || config.excludedGenes.contains(interaction.to.name))))

    println("Interactions: " + interactions.size)

    println("Loading PValues file...")
    // Read lines while dropping the header
    val pValuesMap = Source.fromFile(config.input_list).getLines.drop(1).map(line => {
      val splitted = line.split(",")
      (splitted(0), splitted(1).toDouble)
    }).toMap

    // Select genes to be used as seed genes
    val geneList = new HashSet[Gene]
    network.genes.foreach(gene => {
      if (pValuesMap.getOrElse(gene.name, 0d) > config.seedGenesMutations) {
        geneList.add(gene)
      }
    })

    (interactions, network, pValuesMap, geneList.view.toSet)
  }
  
}