package be.cmpg.cancer

import be.cmpg.walk.fungus.Fungus
import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.graph.Interaction

class Hi2iiFungus(startGene: Gene,
  endGenes: Set[Gene] = Set(),
  geneNumberVariable: Double = 3,
  network: NetworkManager[_]) extends Fungus(startGene, endGenes, geneNumberVariable, network) {

}