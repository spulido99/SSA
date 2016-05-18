package be.cmpg.graph.interaction

import scala.collection.Set
import scala.collection.Map
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.Gene
import be.cmpg.walk.Path
import scala.collection.Iterable
import scala.collection.Traversable
import be.cmpg.walk.SubNetworkSelector
import scala.util.control.Breaks
import scala.util.Random
import be.cmpg.graph.Interaction

/*
 * The cost belong to the interactions.
 */
class InteractionCostNetworkManager(
  network: Network,
  pheromone: Double = 1.0005,
  evaporation: Double = 0.999) extends NetworkManager[Interaction](network) {

  val probabilityMap = initProbabilityMap()

  def initProbabilityMap() = {
    val map = new scala.collection.mutable.HashMap[Interaction, Double]()
    network.interactions.map(interaction => map(interaction) = 0.5d)
    map
  }

  /*
   * a subnetowrk does not have a specific default score, its score is given by the number of times it is used
   */
  def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None) : Double = Double.NaN
  
  override def updateScores(subnetworkScores: Traversable[WalkerResult]) : Map[Interaction, Double] = {
    subnetworkScores.map(_.subnetwork).flatten.toSet[Interaction].foreach(interaction => probabilityMap.put(interaction, math.min(1, probabilityMap(interaction) * pheromone)))
    probabilityMap.view.toMap
  }

  override def evaporate(geneScores:Map[Interaction, Double]) = {
    probabilityMap.foreach(x => probabilityMap.put(x._1, math.max(0, x._2 * evaporation)))
  }

  override def getPosteriorProbability(interaction: Interaction): Double =
    probabilityMap(interaction) * interaction.probability
    
    /**
   * Selects a random interaction from a list of Interactions based on the posterior probability as defined
   * in the NetworkManager. Thus an Interaction with a high posterior probability will be picked more than one with a 
   * low probability.
   * 
   * @returns None if no interaction can be selected, Some(Interaction) where the Interaction is randomly picked.
   */
  override def getRandomInteraction(selector: SubNetworkSelector): Option[Interaction] = {
   
    val childInteractions = selector.getPossibleInteractions()
    
    if (childInteractions.isEmpty) {
      None
    } else {
      var ran = Random.nextDouble * childInteractions.map(child => getPosteriorProbability(child)).sum
      var selectedInteraction: Interaction = null

      Breaks.breakable {
    	  for (interaction <- childInteractions) {
    	    selectedInteraction = interaction
    	    ran -= getPosteriorProbability(selectedInteraction)
    	    
    	    if (ran <= 0) Breaks.break
    	  }
      } 
      Some(selectedInteraction)
    }
  }

}