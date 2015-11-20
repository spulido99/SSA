package be.cmpg.walk.fungus

import org.junit.runner.RunWith
import org.specs2.mutable.Specification
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import org.specs2.runner.JUnitRunner
import be.cmpg.graph.Network
import be.cmpg.graph.interaction.InteractionCostNetworkManager
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.expression.ExpressionNetworkManager

@RunWith(classOf[JUnitRunner])
class FungusSpecification extends Specification {

  "A Fungus" should {
    "represent a growing fungus, one interaction at a time" in {

      val network = new ExpressionNetworkManager( new Network(
        interactions = Set(
          Interaction(Gene("from"), Gene("to"), probability = 1),
          Interaction(Gene("to"), Gene("end_1"), probability = 1),
          Interaction(Gene("to"), Gene("end_2"), probability = 1))))

      "it should contain only the start node on initialization" in {
        val fungus = new Fungus(Gene("from"), Set(), 1, network)
        fungus.visitedNodes.contains(network.getNetwork().getNode(Gene("from"))) must beTrue
      }

      "it should be expandable if a max value of 2  (genes), then just one interation should be there" in {
        val fungus = new Fungus(Gene("from"), Set(), 2, network)
        val interaction = Interaction(from = Gene("from"), to = Gene("to"))
        
        val path = fungus.selectSubNetwork().get
        path.size must be_== (1)
      }

      "it should finish if the complete fungus overgrows the complete network" in {
        val fungus = new Fungus(Gene("from"), Set(), 1000, network)
        val path = fungus.selectSubNetwork().get
        path.map(_.genes).flatten.size == 4 must beTrue
      }

    }
  }
}