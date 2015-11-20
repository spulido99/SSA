package be.cmog.benchmark

import be.cmpg.graph.NetworkReader
import java.nio.file.Paths
import be.cmpg.graph.interaction.NetworkManager
import scala.io.Source
import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.InteractionCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.walk.fungus.DFungus
import be.cmpg.optimization.AntLikeOptimization

object BenchmarkBasedOnImidazole {

  def generateDFungus(network: NetworkManager[_], genes: Set[Gene]) = {
    (0 to 1).map(_ => genes.map(gene => {
      new DFungus(
        startGene = gene,
        endGenes = genes - gene,
        maximumIterations = 100,

        network = network)
    })).flatten
  }

  def main(args: Array[String]) {

    val interactions = NetworkReader
    		.fromFile(Paths.get("src/main/resources/be/cmpg/benchmark/salmonella_only_metabolic_network.txt"))

    val networkManager = new InteractionCostNetworkManager(
      network = new Network(interactions),
      evaporation = 0.999999d)

    val inputGeneFile = Paths.get("src/main/resources/be/cmpg/benchmark/input_genes_imi_log2_fc_lt_0_75.txt")

    val inputGenes = Source.fromFile(inputGeneFile.toFile)
      .getLines
      .filter(!_.isEmpty())
      .map(Gene(_))
      .toSet
      .filter(networkManager.getGenes.contains(_))
      
    val fungi = generateDFungus(networkManager, inputGenes)
    
    val optimization = new AntLikeOptimization(networkManager, fungi)

    optimization.optimize(100)

    networkManager.printResults(System.out)
  }
}