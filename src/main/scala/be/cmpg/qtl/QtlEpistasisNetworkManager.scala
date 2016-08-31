package be.cmpg.qtl

import org.apache.commons.math3.distribution.BetaDistribution
import be.cmpg.graph.Gene
import be.cmpg.graph.Network
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Interaction
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import scala.collection.Set
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.utils.weightByFlatInitialProbability

class QtlEpistasisNetworkManager(network: Network,
    chromosomeBackground: Map[String, BetaDistribution],
    markerFrequency: Map[Set[String], Double], // tuple: Chr:Pos ChrXI:149594
    pheromone: Double,
    initialProb:Double= 0.5,
    evaporation: Double ) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double,weightingScheme=new weightByFlatInitialProbability(network,initialProb)) {

  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    
    val startLociMarkers = subnetwork.filter(_.typ == QtlConstants.startLoci).map(_.genes).flatten
    val endLociMarkers = subnetwork.filter(_.typ == QtlConstants.endLoci).map(_.genes).flatten
    
    val startLociStats = new SummaryStatistics
    
    for (marker <- startLociMarkers) {
      
      val chr = marker.get(QtlConstants.chromosome).get
      val pos = marker.get(QtlConstants.position).get
      
      startLociStats.addValue(chromosomeBackground(chr).cumulativeProbability(markerFrequency(Set(chr, pos))))
    }
    
    val endLociStats = new SummaryStatistics
    for (marker <- endLociMarkers) {
      
      val chr = marker.get(QtlConstants.chromosome).get
      val pos = marker.get(QtlConstants.position).get
      
      endLociStats.addValue(chromosomeBackground(chr).cumulativeProbability(markerFrequency(Set(chr, pos))))
    }
    
    val startPVal = 0
    val endPVal = 0
    
    /*
     * 
     * How do we account for different distances between markers?
     * 
     * - The size of each loci should be small, so the number of 
     * recombination events should be asumed to be close to zero in the region.
     * There is no need to account for the different distance between markers in each loci
     * 
     * How do we account for different amount of markers in each loci? Average? Median?
     * 
     * Average have the problem to be affected a lot by just one very high frequency marker (should NOT use average)
     * 
     * 
     * 
     * 
     * 
     */ 
    
    math.min(startPVal, endPVal).toDouble
  }
  
}