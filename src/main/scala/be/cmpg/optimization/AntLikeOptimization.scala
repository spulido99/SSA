package be.cmpg.optimization

import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.graph.Interaction
import java.util.concurrent.Callable
import scala.collection.Set
import be.cmpg.graph.interaction.WalkerResult

class AntLikeOptimization(networkManager: NetworkManager[_], subNetworkSelectors: Traversable[SubNetworkSelector]) {

  def optimize(numberOfSteps: Int) {
    val threadpool = java.util.concurrent.Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors())

    for (i <- 0 to numberOfSteps) {

      val callables = subNetworkSelectors.map(ant => {
        new Callable[Option[WalkerResult]] {
          override def call(): Option[WalkerResult] = {
            val subnetwork = ant.selectSubNetwork()
            if (subnetwork.isDefined)
              Some(WalkerResult(ant, subnetwork.get, networkManager.scoreSubnetwork(subnetwork.get)))
            else
              return None
          }
        }
      })

      val previous = System.nanoTime()
      
      val futures = callables.map(threadpool.submit(_))

      val after = System.nanoTime()

      val selectedSubNetwork = futures.map(_.get()).flatten.toSet 

      networkManager.updateScores(selectedSubNetwork)

      if (i != numberOfSteps)
        networkManager.evaporate(Map())

      if (i % 100 == 0) {
        println("Finished pathfinding for round " + i + " in " + (after.toDouble - previous.toDouble) / 1000d + " ms.")
      }

    }

    threadpool.shutdown()
  }

}