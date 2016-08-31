package be.cmpg.walk

import org.junit.runner.RunWith
import org.specs2.mutable.Specification
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import org.specs2.runner.JUnitRunner

@RunWith(classOf[JUnitRunner])
class PathSpecification extends Specification {

  "A Path" should {
    "represent a randomwalk state in the current network" in {

      "it should be empty on initialization" in {
        val path = new Path(startGene = Gene("from"))
        path.getVisitedInteractions() must have size(0)
        path.currentEndpoint must beEqualTo(Gene("from"))
        path.size() must beEqualTo(0)
      }

      "it should be expandable only when it can take an interaction" in {
        val path = new Path(Gene("from"))

        val firstInteraction = Interaction(Gene("from"), Gene("to"), "typ")

        path.canTakeInteraction(firstInteraction) must beTrue
        path.expand(firstInteraction)

        path.getVisitedInteractions must contain(allOf(firstInteraction))
        path.currentEndpoint must beEqualTo(Gene("to"))
        path.size() must beEqualTo(1)
      }

      "it should throw an IllegalArgumentException when an element is added that cannot used to expand the RandomWalk" in {
        val walk = new Path(Gene("from"))

        val firstInteraction = Interaction(Gene("from"), Gene("to"), "typ")
        walk.expand(firstInteraction)

        "when it revisits a node already visited by the RandomWalk" in {
          val secondInteraction = Interaction(Gene("to"), Gene("from"), "typ")

          walk.canTakeInteraction(secondInteraction) must beFalse
          walk.expand(secondInteraction) must throwA[IllegalArgumentException]
        }

        "when an interaction is added which is not connected with the current endpoint" in {
          val secondInteraction = Interaction(Gene("to_1"), Gene("from_1"), "typ")

          walk.canTakeInteraction(secondInteraction) must beFalse
          walk.expand(secondInteraction) must throwA[IllegalArgumentException]
        }
      }

      "it should be resetable" in {
        val walk = new Path(Gene("from"))

        val firstInteraction = Interaction(Gene("from"), Gene("to"), "typ")
        walk.expand(firstInteraction)

        walk.reset

        walk.currentEndpoint must beEqualTo(Gene("from"))
        walk.getVisitedInteractions must have size(0)
        walk.visitedGenes must contain (allOf(Gene("from")))
      }
    }
  }
}