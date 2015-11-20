package be.cmpg.graph

import java.nio.file.Paths
import java.util.concurrent.Callable
import scala.collection.JavaConversions._
import scala.io.Source
import be.cmpg.walk.Ant
import be.cmpg.walk.Path
import be.cmpg.walk.fungus.Fungus
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.graph.interaction.InteractionCostNetworkManager
import be.cmpg.optimization.AntLikeOptimization
import be.cmpg.walk.fungus.DFungus

object Phac {

  val genesAsString = """STM2540
STM2380
STM1173
STM0366
STM2279
STM2771
STM3422
STM1977
STM3648
STM1153
STM1975
STM2782
STM3443
STM1851
STM2801
STM4295
STM2024
STM0543
STM1202
STM0465
STM3033
STM2831
STM1976
STM1803
STM3034
STM4538
STM1959
STM3245
STM4267
STM0218
STM3798
STM4394
STM1564
STM1973
STM1988.S
STM0474
STM2022
STM3687
STM2299
STM1285
STM3414
STM3229
STM3225
STM1831
STM1960
STM4467
STM1246
STM4459
STM2027
STM3630
STM2675
STM1311
STM1898
STM1572
STM3270
STM0087
STM2065
STM1118
STM1832
STM1308
STM2929
STM2023
STM4240
STM2311
STM4460
STM1731
STM3471
STM0141
STM3384
STM0421
STM0687
STM3325
STM2673
STM4465
STM3385
STM1978
STM2031
STM0701
STM4519
STM1563
STM4561
STM4026
STM2796
STM0665
STM1119
STM2021
STM3435
STM2369
STM2297
STM1318
STM2026
STM1255
STM1121
STM1984
STM4229
STM2802
STM4392
STM1223
STM1796
STM3335
STM0389
STM2141
STM0211
STM0549
STM1376
STM4377
STM0756
STM4511
STM3344"""

  val interactions = NetworkReader.fromFile(Paths.get("src/test/resources/be/cmpg/graph/network_1.txt"))

  def generateAnts(network: NetworkManager[Interaction]) = {
    val genes = {
      Source.fromString(genesAsString).getLines.map(geneName => Gene(geneName.toLowerCase)).toSet.filter(network.getGenes.contains(_))
    }
    (0 to 0).map(_ => genes.map(gene => {
      new Ant(
        restartFrequency = 0.15d,
        maximumIterations = 5000,

        startGene = gene,
        endGenes = genes - gene,

        network = network)
    })).flatten
  }

  def generateFungus(network: NetworkManager[Gene], score: Double) = {
    val genes = {
      Source.fromString(genesAsString).getLines.map(geneName => Gene(geneName.toLowerCase)).toSet.filter(network.getGenes.contains(_))
    }
    (0 to 0).map(_ => genes.map(gene => {
      new Fungus(
        startGene = gene,
        endGenes = genes - gene,

        _geneNumberVariable = score,

        network = network)
    })).flatten
  }

  def generateDFungus(network: NetworkManager[_], genes: Set[Gene]) = {
    (0 to 1).map(_ => genes.map(gene => {
      new DFungus(
        startGene = gene,
        endGenes = genes - gene,
        maximumIterations = 100,

        network = network)
    })).flatten
  }

  def loadInputGenes(network: NetworkManager[_]) = {
    Source.fromString(genesAsString)
      .getLines.map(geneName => Gene(geneName.toLowerCase))
      .toSet
      .filter(network.getGenes.contains(_))
  }

  def main(args: Array[String]) {

    val networkManager = new InteractionCostNetworkManager(
      network = new Network(interactions),
      evaporation = 0.999999d)

    val inputGenes = loadInputGenes(networkManager)

    val steps = 1000
    val ants = generateAnts(networkManager)
    val dfungi = generateDFungus(networkManager, inputGenes)
    //val walkers = generateFungus(network, 20)

    val optimization = new AntLikeOptimization(networkManager, dfungi)

    optimization.optimize(steps)

    networkManager.printResults(System.out)

  }

}