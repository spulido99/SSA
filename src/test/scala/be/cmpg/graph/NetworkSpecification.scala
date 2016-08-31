package be.cmpg.graph

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import java.nio.file.Paths

@RunWith(classOf[JUnitRunner])
class NetworkSpecification extends Specification {

  "A Network" should {

    "be able to find the outgoing nodes of a given gene" in {
      val network = new Network(Set(
        Interaction(Gene("from"), Gene("to"), direction = "undirected"),
        Interaction(Gene("to"), Gene("end_1"), direction = "undirected"),
        Interaction(Gene("to"), Gene("end_2"), direction = "undirected")))
      

      network.getOutgoingInteractions(Gene("to")) must contain(allOf(
        Interaction(Gene("to"), Gene("from")),
        Interaction(Gene("to"), Gene("end_1")),
        Interaction(Gene("to"), Gene("end_2"))))
    }
  }
}