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
import be.cmpg.utils.testWeighting
import be.cmpg.utils.weightByFlatInitialProbability

@RunWith(classOf[JUnitRunner])
class FungusSpecification extends Specification {

  "A Fungus" should {
    "represent a growing fungus, one interaction at a time" in {

      val networkFlatWeighted = new Network(
        interactions = Set(
          Interaction(Gene("from"), Gene("to")),
          Interaction(Gene("to"), Gene("end_1")),
          Interaction(Gene("to"), Gene("end_2")),
          Interaction(Gene("to"), Gene("end_3"))))

      val networkManager = new ExpressionNetworkManager(networkFlatWeighted, weightingScheme = new weightByFlatInitialProbability(networkFlatWeighted, 0.5))

      val networkEnd_2Downweighted = new Network(
        interactions = Set(
          Interaction(Gene("from"), Gene("to")),
          Interaction(Gene("to"), Gene("end_1")),
          Interaction(Gene("to"), Gene("end_2")),
          Interaction(Gene("to"), Gene("end_3"))))

      val networkEnd_3Downweighted = new Network(
        interactions = Set(
          Interaction(Gene("from"), Gene("to")),
          Interaction(Gene("to"), Gene("end_1")),
          Interaction(Gene("to"), Gene("end_2")),
          Interaction(Gene("to"), Gene("end_3")),
          Interaction(Gene("end_1"), Gene("end_5")),
          Interaction(Gene("end_1"), Gene("end_4"))))

      val downweightedNetworkManager_3 = new ExpressionNetworkManager(networkEnd_3Downweighted, weightingScheme = new testWeighting(networkEnd_3Downweighted, 0.5, List("end_1"), 0))
      
      val downweightedNetworkManager = new ExpressionNetworkManager(networkEnd_2Downweighted, weightingScheme = new testWeighting(networkEnd_2Downweighted, 0.5, List("end_2"), 0))

      "it should contain the start node when having selected a subnetwork" in {
        val fungus = new Fungus(Gene("from"), Set(), 1, networkManager)
        fungus.selectSubNetwork()
        fungus.visitedNodes.contains(networkFlatWeighted.getNode(Gene("from"))) must beTrue
      }

      "it should be expandable if a max value of 2  (genes), then just one interation should be there" in {
        val fungus = new Fungus(Gene("from"), Set(), 2, networkManager)
        val interaction = Interaction(from = Gene("from"), to = Gene("to"))
        val path = fungus.selectSubNetwork().get
        path.size must be_==(1)
      }

      "it should finish if the complete fungus overgrows the complete network" in {
        val fungus = new Fungus(Gene("from"), Set(), 1000, networkManager)
        val path = fungus.selectSubNetwork().get
        path.map(_.genes).flatten.size == 5 must beTrue
      }

      "it should 'randomly' select genes based on their initial probability. So an interaction with 0 initial probability should never be selected" in {
        val fungus = new Fungus(Gene("from"), Set(), 100, downweightedNetworkManager)
        val selectedInteractions = fungus.selectSubNetwork().get
        val selectedGenes = selectedInteractions.map(interaction => interaction.genes).flatten.map(gene => gene.name)
        selectedGenes should not contain ("end_2")
      }

      "it should select a subnetwork which is defined as a connected set of edges" in {
        val fungus_1 = new Fungus(Gene("from"), Set(), 100, downweightedNetworkManager_3)
        val selectedSubnetwork = fungus_1.selectSubNetwork().get
        var connected = true
        selectedSubnetwork.foreach(interaction => {
          val otherInteractions = selectedSubnetwork.-(interaction)
          val otherGenes = otherInteractions.toList.map(interaction => interaction.genes).flatten.map(gene => gene.name)
          if (!(otherGenes.contains(interaction.from.name)) && !(otherGenes.contains(interaction.to.name))){connected = false}
        })
        connected should beEqualTo(true)
      }
      
    }
  }
}