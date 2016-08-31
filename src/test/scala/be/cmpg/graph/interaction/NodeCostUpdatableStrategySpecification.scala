package be.cmpg.graph.interaction

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import be.cmpg.graph.Gene
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.walk.Ant
import be.cmpg.walk.Path
import be.cmpg.walk.fungus.SubNetwork
import be.cmpg.expression.ExpressionNetworkManager
import be.cmpg.utils.weightByFlatInitialProbability

@RunWith(classOf[JUnitRunner])
class NodeCostUpdatableStrategySpecification extends Specification {

  "An NodeCostUpdatableStrategy" should {
    "be able to increase the value of used Genes and decrease the value of Unused genes" in {

      val network = new Network(
          interactions = Set(
            Interaction(Gene("from"), Gene("to")),
            Interaction(Gene("from"), Gene("useless_1")),
            Interaction(Gene("from"), Gene("useless_2")),
            Interaction(Gene("useless_1"), Gene("useless_3")),
            Interaction(Gene("useless_2"), Gene("useless_3")),
            Interaction(Gene("to"), Gene("end_1")),
            Interaction(Gene("to"), Gene("end_2"))))
      
      // take evaporation 1 so nothing evaporates
      val interactionManager = new ExpressionNetworkManager(
       network, evaporation=1.0,weightingScheme = new weightByFlatInitialProbability(network,0.5))

      val path = new SubNetwork(Gene("from"))
      path.expand(Interaction(Gene("from"), Gene("to")))
      path.expand(Interaction(Gene("to"), Gene("end_1")))
      
      // Specific values are irrelevant as the subnetwork is imposed beforehand
      val expressionData = Map((Gene("from"),0.5),(Gene("to"),0.5),(Gene("useless_1"),0.5),(Gene("useless_2"),0.5),(Gene("useless_3"),0.5),(Gene("end_1"),0.5),(Gene("end_2"),0.5))
      
      interactionManager.geneExpression=expressionData
      
      val subnetwork = path.getVisitedInteractions().toSet
      val score = interactionManager.scoreSubnetwork(subnetwork)
      
      interactionManager.updateScores(List((subnetwork, score)))
      interactionManager.evaporate()

      // Check that the current value of "from" and "end_x" -the found end- node are higher than the value of "to"
      // Check that the value of useless_n are lower than the original 0.5
      // (it can be that they are part of the path, righ?)
      "the score from the genes in the path should be higher than 0.5" in {
        interactionManager.getPosteriorProbability(Gene("to")) must be_> (0.5)
      }
    }
  }
}