package be.cmpg.walk.neighbours

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import be.cmpg.graph.Gene
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.InteractionCostNetworkManager

@RunWith(classOf[JUnitRunner])
class NeighbourhoodWatchSpecification extends Specification {

  "A NeighbourhoodWatch" should {
    "be able to find the interactions which are located near a given gene" in {

      val startGene = new Gene("from")

      val watch = new NeighbourhoodWatch(startGene = startGene,
        depth = 1,
        network = new InteractionCostNetworkManager(new Network(
          interactions = Set(
            Interaction(Gene("from"), Gene("to"), probability = 1),
            Interaction(Gene("end_1"), Gene("end_3"), probability = 1),
            Interaction(Gene("to"), Gene("end_1"), probability = 1),
            Interaction(Gene("to"), Gene("end_2"), probability = 1)))))

      val result = watch.selectSubNetwork()

      result must beSome

      result.get.size must beEqualTo(3)
    }
  }
}