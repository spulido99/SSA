package be.cmpg.walk

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import be.cmpg.graph.Gene
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.InteractionCostNetworkManager

@RunWith(classOf[JUnitRunner])
class AntSpecification extends Specification {

  "An Ant" should {
    "be able to walk the interaction network" in {
      val ant = new Ant(
        restartFrequency = 0.33,
        maximumIterations = 100,
        startGene = Gene("from"),
        endGenes = Set(Gene("end_1"), Gene("end_2")),
        network = new InteractionCostNetworkManager(new Network(
          interactions = Set(
            Interaction(Gene("from"), Gene("to")),
            Interaction(Gene("to"), Gene("end_1")),
            Interaction(Gene("to"), Gene("end_2"))))))

      val path = ant.selectSubNetwork()

      "the resulting path should connect the start gene to one of the end genes" in {
        path.get.map(_.genes).flatten.intersect(Set(Gene("end_1"), Gene("end_2"))) must haveSize(1)
      }
    }
  }
}