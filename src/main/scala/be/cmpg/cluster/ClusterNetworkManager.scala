package be.cmpg.cluster

import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import scala.collection.Map
import scala.collection.Set
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.utils.weightByDegreeExponentialdistribution
import be.cmpg.utils.weightByFlatInitialProbability


class ClusterNetworkManager(network: Network,
  genePValueMatrix: Map[String, Double], // (EntrezGeneId, PatientId) -> chrom:start-end
  pheromone: Double = 0.005,
  evaporation: Double = 0.996,
  initialProb:Double = 0.5,
  ranked: Boolean = false,
  minProb: Double = 0.01) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, new weightByDegreeExponentialdistribution(network,initialProb), ranked: Boolean, minProb: Double, initialProb= initialProb) {
  
  println("MEManager")
  println("Genes: "+genePValueMatrix.size)
  
  override def scoreWalker(walker: SubNetworkSelector) : Option[(Set[Interaction], Double)] = {
    var chances = 5
    var result = (Set[Interaction](), Double.NaN)
    while (chances > 0 && result._2.isNaN()) {
      val subnetwork = walker.selectSubNetwork()
      result = (subnetwork.get, scoreSubnetwork(subnetwork.get, Some(walker)))
      chances -= 1
    }
    Some(result)
  }
  
  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    val allGenes = subnetwork.map(_.genes).flatten.toSet

     // Check at least 2 genes are present in the subnetwork
     
    if (allGenes.size < 2)
      return Double.NaN;
      
       
      // Make a list of samples to be analysed (samples with mutations in those genes)
      
      val scorePerGene = allGenes.map(gene => genePValueMatrix(gene.name))

       return scorePerGene.sum / allGenes.size
    }
  
}
 