package be.cmpg.graph
import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import java.nio.file.Paths
import be.cmpg.cluster.NetworkFactory
import java.io.File

@RunWith(classOf[JUnitRunner])
class NetworkDegreeSpecification extends Specification {
  
  
  "the cluster NetworkFactory" should{
    "be able to read a network in simple csv format" in{
        val network = NetworkFactory.loadNetwork(new File("src/main/resources/networks/clustering/net.csv"))
        network.interactions.size must beGreaterThan(1)
    }
  }
}