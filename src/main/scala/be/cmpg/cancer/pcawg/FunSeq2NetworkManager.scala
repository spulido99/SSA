package be.cmpg.cancer.pcawg

import be.cmpg.graph.Network
import scala.collection.Set
import be.cmpg.graph.Gene
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Interaction
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.walk.fungus.Fungus
import scala.util.Try

class FunSeq2NetworkManager(network: Network,
                              geneFunSeqScore: List[FunSeqData],
                              pheromone: Double,
                              evaporation: Double,
                              ranked: Boolean,
                              convergenceThreshold:Double) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, ranked: Boolean, convergenceThreshold:Double) {

	val bySample = geneFunSeqScore.groupBy { _.sample }.mapValues { _.groupBy { _.gene }.mapValues { _.maxBy { _.score } } };
  val hypermutators = bySample.filter(e => e._2.size > 400).keySet
  val byGene   = geneFunSeqScore.groupBy { _.gene }.mapValues { _.groupBy { _.sample }.mapValues { _.maxBy { _.score } }.toList.filterNot(e => hypermutators.contains(e._1)) };
  //val byGene   = geneFunSeqScore.groupBy { _.gene }.mapValues { _.groupBy { _.sample }.mapValues { _.maxBy { _.score } }.toList };
  
  {
	  println("Samples: "+bySample.size)
	  println("Genes : "+ byGene.size)
  }
  
  
  def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    val genes = subnetwork.map(_.genes).flatten
    val startGene = if (selector.isDefined) selector.get.getStartGene() else Gene("")
    
    
    
    if (selector.isDefined && genes.size < selector.get.getGeneNumberVariable) return Double.NaN
    
    /*
     * The score of a subnetwork is equal to the min of the 
     */
    
    val geneScores = genes.map { g => byGene.getOrElse(g, List()).foldLeft(0.0)(_ + _._2.score) }
    
    geneScores.sum
    
    /*
     * the max score per sample in the genes
     *
    val scores = bySample.values.map { mutationsBySample => 
         genes.map { g =>
           mutationsBySample.get(g).map { _.score }
         }
       }.flatten.flatten
       
    val meanScore = if (scores.isEmpty) 0.0 else scores.reduce(_ + _)/scores.size
       //.maxBy { _.score }
       //.score
    
    /
    var result = 0.0
    var resultMinusStart = 0.0
    genes.foreach{ gene => {
         val score = geneFunSeqScore.get(gene)
         result += score.getOrElse(0.0)
         if (gene != startGene) {
           resultMinusStart += score.getOrElse(0.0)
         }
      }}
      
    result /= genes.size
    resultMinusStart /= genes.size - 1
    // Return NaN if the network would be better off without the starting gene
    if (resultMinusStart < result) Double.NaN else result*/
    
    //meanScore
  }

}