package be.cmpg.walk.fungus

import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.Gene
import java.util.concurrent.Callable
import scala.collection.Set
import be.cmpg.expression.ExpressionNetworkManager
import be.cmpg.utils.weightByFlatInitialProbability

object FungusTestFixme {

  val interactions = Set(
    Interaction(Gene("1"), Gene("2")),
    Interaction(Gene("2"), Gene("3")),
    Interaction(Gene("3"), Gene("4")),
    Interaction(Gene("4"), Gene("5")),
    Interaction(Gene("4"), Gene("6")),
    Interaction(Gene("4"), Gene("7")),
    Interaction(Gene("6"), Gene("8")),
    Interaction(Gene("7"), Gene("8")),
    Interaction(Gene("7"), Gene("10")),
    Interaction(Gene("8"), Gene("9")),
    Interaction(Gene("10"), Gene("11")),
    Interaction(Gene("11"), Gene("12")),
    Interaction(Gene("11"), Gene("14")),
    Interaction(Gene("12"), Gene("14")),
    Interaction(Gene("12"), Gene("13")),
    Interaction(Gene("13"), Gene("18")),
    Interaction(Gene("14"), Gene("15")),
    Interaction(Gene("15"), Gene("16")),
    Interaction(Gene("16"), Gene("17")),
    Interaction(Gene("17"), Gene("18")))
  val network = new Network(interactions)
  val networkManager = new ExpressionNetworkManager(network = network, weightingScheme= new weightByFlatInitialProbability(network,0.5))
  
  val notWantedScore = 10;
  val wantedScore = 3;
  
  network.getNode(Gene("1")).score = notWantedScore
  network.getNode(Gene("2")).score = notWantedScore
  network.getNode(Gene("3")).score = wantedScore
  network.getNode(Gene("4")).score = wantedScore
  network.getNode(Gene("5")).score = notWantedScore
  network.getNode(Gene("6")).score = wantedScore
  network.getNode(Gene("7")).score = notWantedScore
  network.getNode(Gene("8")).score = wantedScore
  network.getNode(Gene("9")).score = wantedScore
  network.getNode(Gene("10")).score = notWantedScore
  network.getNode(Gene("11")).score = notWantedScore
  network.getNode(Gene("12")).score = wantedScore
  network.getNode(Gene("13")).score = wantedScore
  network.getNode(Gene("14")).score = wantedScore
  network.getNode(Gene("15")).score = wantedScore
  network.getNode(Gene("16")).score = notWantedScore
  network.getNode(Gene("17")).score = notWantedScore
  network.getNode(Gene("18")).score = notWantedScore


  def generateFungus(network: NetworkManager[Gene], score: Double) = {
    (0 to 0).map(_ => networkManager.getGenes.map(gene => {
      new Fungus(
        startGene = gene,
        endGenes = Set(),

        _geneNumberVariable = score,

        network = network)
    })).flatten
  }

  def main(args: Array[String]) {
    //val threadpool = java.util.concurrent.Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors())
    val threadpool = java.util.concurrent.Executors.newSingleThreadExecutor()

    val numberOfSteps = 1000

    for (i <- 0 to numberOfSteps) {
      
      val callables = generateFungus(networkManager, 19).map(walker => {
        new Callable[Option[(Set[Interaction], Double)]] {
          override def call(): Option[(Set[Interaction], Double)] = {
            val subnetwork = walker.selectSubNetwork()
            if (subnetwork.isDefined)
              Some((subnetwork.get, networkManager.scoreSubnetwork(subnetwork.get)))
            else
              return None
          }
        }
      })

      val previous = System.nanoTime()

      val futures = callables.map(threadpool.submit(_))

      val after = System.nanoTime()

      //futures.map(x => println(x.get))
      val paths = futures.map(_.get.get)

      networkManager.updateScores(paths)

      if (i != numberOfSteps)
        networkManager.evaporate()

      if (i % 100 == 0) {
        println("Finished pathfinding for round: " + i + " in: " + (after.toDouble - previous.toDouble) / 1000d + " ms.")
      }

    }

    threadpool.shutdown()

    networkManager.printResults(System.out)

  }

}