package be.cmpg.walk

import be.cmpg.graph.Gene
import be.cmpg.walk.fungus.Hi2iiFungus
import be.cmpg.utils.StatUtils
import be.cmpg.graph.interaction.NodeCostNetworkManager

object WalkerFactory {
  
def buildHi2iiFungusWalker(geneList: Set[Gene], networkManager: NodeCostNetworkManager) = {
    val walkers: Set[SubNetworkSelector] = geneList.map(gene =>
      //Initialize new Fungus
      new Hi2iiFungus(
        startGene = gene,
        network = networkManager)).toSet

    walkers
  }
  
}