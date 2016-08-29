package be.cmpg.walk.fungus

import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.NetworkManager

class Hi2iiFungus(startGene: Gene,
  endGenes: Set[Gene] = Set(),
  geneNumberVariable: Double = 0,
  network: NetworkManager[_]) extends Fungus(startGene, endGenes, geneNumberVariable, network) {

}