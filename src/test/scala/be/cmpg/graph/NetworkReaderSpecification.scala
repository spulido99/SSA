package be.cmpg.graph

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import java.nio.file.Paths

@RunWith(classOf[JUnitRunner])
class NetworkReaderSpecification extends Specification{

  "The GraphReader" should {
    "convert a given text file to a graph" in {
      
      val graph = NetworkReader.fromFile("src/test/resources/be/cmpg/graph/network_1.txt")
    
      graph.size must beEqualTo(13327)
      
    }
  }
}