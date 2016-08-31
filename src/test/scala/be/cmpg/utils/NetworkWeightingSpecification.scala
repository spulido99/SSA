package be.cmpg.utils

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import be.cmpg.cluster.NetworkFactory
import java.io.File

@RunWith(classOf[JUnitRunner])
class NetworkWeightingSpecification extends Specification {
  
    "A Netzorkzeighting" should {
      "Be able to initialize node probabilities on the network" in{
  
  val network = NetworkFactory.loadNetwork(new File("src/main/resources/networks/clustering/net.csv"))
  val weightedNetworkByFlatDistribution = new weightByFlatInitialProbability(network,0.5).initialize()
  
  val probabilityMap = weightedNetworkByFlatDistribution.nodeMap.map(node => (node._1.name,node._2.posteriorProbability))
  probabilityMap("2064") should beEqualTo(0.5)
    }
    }
}