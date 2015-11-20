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
        Interaction(Gene("from"), Gene("to"), direction = "undirected", probability = 1d),
        Interaction(Gene("to"), Gene("end_1"), direction = "undirected", probability = 1d),
        Interaction(Gene("to"), Gene("end_2"), direction = "undirected", probability = 1d)))

      network.getOutgoingInteractions(Gene("to")) must contain(allOf(
        Interaction(Gene("from"), Gene("to"), probability = 1d),
        Interaction(Gene("to"), Gene("end_1"), probability = 1d),
        Interaction(Gene("to"), Gene("end_2"), probability = 1d)))
    }
  }
}