package be.cmpg.walk.fungus

import org.junit.runner.RunWith
import org.specs2.mutable.Specification
import org.specs2.runner.JUnitRunner
import be.cmpg.graph.Gene
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.InteractionCostNetworkManager

@RunWith(classOf[JUnitRunner])
class DFungusSpecification extends Specification {

  "A DFungus" should {
    "select a set of interactions of the network which is in essence a search tree that starts in the start gene and ends in target genes" in {
      val startGene = Gene("start")

      val endGene_1 = Gene("end-1")
      val endGene_2 = Gene("end-2")

      val endGenes = Set(endGene_1, endGene_2)

      val fungus = new DFungus(
        maximumIterations = 100,
        startGene = Gene("from"),
        endGenes = Set(Gene("end_1"), Gene("end_2")),
        network = new InteractionCostNetworkManager(network = new Network(
          interactions = Set(
            Interaction(Gene("from"), Gene("to"), probability = 1),
            Interaction(Gene("from"), Gene("useless_1"), probability = 1),
            Interaction(Gene("useless_2"), Gene("useless_1"), probability = 1),
            Interaction(Gene("to"), Gene("end_1"), probability = 1),
            Interaction(Gene("to"), Gene("end_2"), probability = 1)))))

      val interactionSet = fungus.selectSubNetwork()

      interactionSet must not beEmpty

      interactionSet.get must contain(exactly(
        Interaction(Gene("from"), Gene("to"), probability = 1),
        Interaction(Gene("to"), Gene("end_1"), probability = 1),
        Interaction(Gene("to"), Gene("end_2"), probability = 1)))
    }
  }
}