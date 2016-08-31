package be.cmpg.expression.time

import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import scala.collection.Set
import be.cmpg.graph.Gene
import com.sun.org.apache.xalan.internal.xsltc.compiler.Expression
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.utils.weightByFlatInitialProbability

class TimeExpressionNM(network: Network,
  timeExpressionMatrix: Map[String, Array[Array[Double]]], // gene -> Matrix[5 Samples][5 Time Points]
  all_samples: Set[String], // all the samples names (here to not need to analyse all keys of the genePatientMatrix)
  mAS_perGene: Int, // minFreqGene: minimum frequency of mutation for a gene to be considered in the analysis
  pheromone: Double = 0.05,
  initialProb:Double = 0.5,
  evaporation: Double = 0.996) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, weightingScheme= new weightByFlatInitialProbability(network,initialProb)) {

  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    var score = 0.0;
    
    if (subnetwork.isEmpty)
      return score
     
    val genes = subnetwork.map(_.genes).flatten[Gene]
    
    genes.foreach(gene => {
      val data = timeExpressionMatrix.get(gene.name)
      if (data.isDefined) {
        score += data.get.flatten.sum
      }
    })

    score/genes.size
  }

}