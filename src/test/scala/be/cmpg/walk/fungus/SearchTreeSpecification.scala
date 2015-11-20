package be.cmpg.walk.fungus

import org.junit.runner.RunWith
import org.specs2.mutable.Specification
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import org.specs2.runner.JUnitRunner

@RunWith(classOf[JUnitRunner])
class SearchTreeSpecification extends Specification {

  "A SearchTree" should {
    "represent a random search from a given point" in {

      "it should contain only the start node on initialization" in {
        val tree = new SearchTree(Gene("start"), Set(Gene("end-1"), Gene("end-2")))
        tree.isGeneVisited(Gene("start"))
      }

      "it should be expandable" in {

        "if an interaction links to a current gene in the search tree" in {
          val tree = new SearchTree(Gene("start"), Set(Gene("end-1"), Gene("end-2")))
          val interaction = Interaction(from = Gene("start"), to = Gene("end"))
          tree.expand(interaction)
          tree.contains(interaction) must beTrue
        }

        "but if an interaction does not link to a current gene in the search tree it should throw an IllegalArgumentException" in {
          val tree = new SearchTree(Gene("start"), Set(Gene("end-1"), Gene("end-2")))
          val interaction = Interaction(from = Gene("some-false-start"), to = Gene("end"))
          tree.expand(interaction) must throwA[IllegalArgumentException]
        }
      }

      "it should be able to return the current leaf genes" in {
        val tree = new SearchTree(Gene("start"), Set(Gene("end-1"), Gene("end-2")))

        tree.expand(Interaction(from = Gene("start"), to = Gene("inter-1")))
        tree.expand(Interaction(from = Gene("inter-1"), to = Gene("end-1")))

        tree.expand(Interaction(from = Gene("inter-1"), to = Gene("end-2")))

        tree.getLeaves must contain(allOf(Set(Gene("end-1"), Gene("end-2"))))
      }

      "it should be trimmable" in {
        val tree = new SearchTree(Gene("start"), Set(Gene("end-1"), Gene("end-2")))

        tree.expand(Interaction(from = Gene("start"), to = Gene("inter-1")))
        tree.expand(Interaction(from = Gene("inter-1"), to = Gene("end-1")))

        tree.expand(Interaction(from = Gene("inter-1"), to = Gene("end-2")))

        tree.expand(Interaction(from = Gene("start"), to = Gene("inter-2")))
        tree.expand(Interaction(from = Gene("inter-2"), to = Gene("inter-3")))
        tree.expand(Interaction(from = Gene("inter-2"), to = Gene("inter-4")))

        val resultingInteractions = tree.trim()

        resultingInteractions must contain(allOf(Set(
          Interaction (from = Gene("start"), to = Gene("inter-1")),
          Interaction(from = Gene("inter-1"), to = Gene("end-2")),
          Interaction(from = Gene("inter-1"), to = Gene("end-1")))))
      }

    }
  }
}