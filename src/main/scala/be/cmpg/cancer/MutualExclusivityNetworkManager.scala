package be.cmpg.cancer

import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import scala.collection.mutable.HashSet
import com.sun.org.apache.xml.internal.serializer.ToSAXHandler
import be.cmpg.graph.Gene
import scala.collection.Map
import scala.collection.Set
import scala.collection.mutable.TreeSet
import java.util.LinkedList
import scala.collection.mutable.HashMap
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import be.cmpg.walk.SubNetworkSelector
import scala.collection.mutable.Buffer
import scala.collection.mutable.MutableList
import scala.collection.mutable.LinearSeq
import scala.collection.immutable.TreeMap
import org.apache.commons.math3.distribution.HypergeometricDistribution
import be.cmpg.graph.interaction.WalkerResult

/**
 * References:
 * An exact algorithm for the weighed mutually exclusive maximum set cover problem
 * http://arxiv.org/abs/1401.6385
 *
 * De novo discovery of mutated driver pathways in cancer
 * http://genome.cshlp.org/content/22/2/375.long
 * The Maximum Coverage Exclusive Submatrix Problem & The Maximum Weight Submatrix Problem
 *
 * Adapted from Dendrix.py from Dendrix_v0.3
 */

class MutualExclusivityNetworkManager(network: Network,
  genePatientMatrix: Map[PolimorphismKey, Polimorphism], // (EntrezGeneId, PatientId) -> chrom:start-end
  minimumSamplesAltered: Int, // minFreqGene: minimum frequency of mutation for a gene to be considered in the analysis
  maximumSamplesAltered: Int = Int.MaxValue,
  pheromone: Double = 0.05,
  evaporation: Double = 0.996,
  ranked: Boolean = false,
  minProb: Double = 0.01,
  statistical: Boolean = false) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, ranked: Boolean, minProb: Double) {
  
  println("MEManager")
  val all_samples = genePatientMatrix.map(_._1.sample).toSet//Set[String], // all the samples names (here to not need to analyse all keys of the genePatientMatrix)
  val mutationsPerSample = genePatientMatrix.groupBy(_._1.sample).map(e => (e._1, e._2.map(e => (e._1.gene, e._2))))
  println("Samples: "+mutationsPerSample.size)
  val mutationsPerGene = genePatientMatrix.groupBy(_._1.gene).map(e => (e._1, e._2.map(e => (e._1.sample, e._2))))
  println("Genes: "+mutationsPerGene.size)
  val mutatedSamplesByGene = mutationsPerGene.keys.map( g => (g, mutationsPerGene(g).keySet)).toMap

  override def scoreWalker(walker: SubNetworkSelector) : Option[WalkerResult] = {
    var chances = 5
    var result = WalkerResult(walker, Set[Interaction](), Double.NaN)
    while (chances > 0 && result.score.isNaN()) {
      val subnetwork = walker.selectSubNetwork()
      result = WalkerResult(walker, subnetwork.get, scoreSubnetwork(subnetwork.get, Some(walker)))
      chances -= 1
    }
    Some(result)
  }
  
  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {
    val allGenes = subnetwork.map(_.genes).flatten.toSet

    val genes = allGenes.filter { g => {
                                	val mutationsInGene = mutationsPerGene.getOrElse(g, Map()).size
                                  (minimumSamplesAltered < mutationsInGene && mutationsInGene < maximumSamplesAltered)
                                }}
    
    /*
     * Check at least 2 genes are present in the subnetwork
     */
    if (genes.size < 2)
      return Double.NaN;

    /*
    if (genes == Set(Gene("MAP2K4"), Gene("MAP3K1"))) {
      avgScores.foreach(s => println(s._1+ " - "+s._2.getMin() + " < "+s._2.getMean() + " < "+s._2.getMax()))
    }
    * 
    */
    
    /*
     * Calculate a score:
     */
    
    // Order genes by to analyse those with more mutations first
    val orderedGenes = genes.toList.sortBy( g => -1*mutatedSamplesByGene(g).size )
    
    // Obtain a score per gene
    /*
     * Score:
     * 
     * Mutation matrix
     *          g1   g2   g3   g4  
     * Sample1   x                       
     * Sample2   x         x    x 
     * Sample3   x                
     * Sample4   x                
     * Sample5   x    x          
     * Sample6        x    x       
     * Sample7        x          
     * Sample8        x           
     * Sample9             x      
     * Sample10                 x
     *      
     */

    if (statistical) {
      /*
       * HyperGeometric test
       */
      //println("********")
      val scorePerGeneHyperTest = orderedGenes.map {
        gene => {
          val otherGenes = genes - gene
          
          val totalSamples = all_samples.size // population size
          val noMutatedSamplesInOtherGenes = all_samples.size - otherGenes.map { mutatedSamplesByGene(_) }.flatten.size // Successes in the population
          val geneMutatedSamples = mutationsPerGene(gene).size //samples size
          
          val test = new HypergeometricDistribution(totalSamples, noMutatedSamplesInOtherGenes, geneMutatedSamples)
          
          val samplesOnlyMutatedInGene = mutationsPerGene(gene).size - mutatedSamplesByGene(gene).intersect(otherGenes.map { mutatedSamplesByGene(_) }.flatten).size // number of succeses
          //println(gene.name + " => PS: "+totalSamples+ " SuccP: "+noMutatedSamplesInOtherGenes+" Sample: "+geneMutatedSamples +" SuccS: " + samplesOnlyMutatedInGene + " => P: "+test.cumulativeProbability(samplesOnlyMutatedInGene))
          
          1.0 - test.cumulativeProbability(samplesOnlyMutatedInGene)
        }
      }
      return scorePerGeneHyperTest.min
      
    } else {
      /*
       * Similar to Dendrix
       */
      // Make a list of samples to be analysed (samples with mutations in those genes)
      val pendingSamples = new HashSet ++ genes.map { g => mutatedSamplesByGene(g) }.reduceLeft(_ ++ _)
      
      val scorePerGene = orderedGenes.map {
        gene => {
          var geneScore = 0.0
          val usedSamples = new HashSet[String]
          for (sample <- pendingSamples) {
            // Mutated samples should be 1 for perfect mutual exclusivity
            // The number samples that have the gene mutated
            if (mutatedSamplesByGene(gene).contains(sample)) {
          	  val mutatedSamples = orderedGenes.map { otherGene => if (mutatedSamplesByGene(otherGene).contains(sample)) 1.0 else 0.0 }.reduceLeft(_ + _)
              usedSamples += sample
              //geneScore += (genes.size - mutatedSamples) / genes.size
              geneScore += 1.0/mutatedSamples
            }
          }
          pendingSamples --= usedSamples
          (gene, geneScore)
        }}
        .toMap
        
       return scorePerGene.map {sg => math.sqrt(sg._2)}.sum / allGenes.size
    }
  }
  
  def cometScore(genes:List[Gene], samples:Set[String]) {
    // not able to replicate
  }
  
  def nonAlteredScore(gene1:Gene, gene2:Gene, nonAlteredSamples:Set[String]):Int = {
    val altG1 = mutatedSamplesByGene(gene1)
    val altG2 = mutatedSamplesByGene(gene2)
    (nonAlteredSamples -- (altG1.diff(altG2) union altG2.diff(altG1))).size
  }
  
  def alteredScore(gene1:Gene, gene2:Gene, alteredSamples:Set[String]):Int = {
	  (alteredSamples -- mutatedSamplesByGene(gene1) -- mutatedSamplesByGene(gene2)).size
  }
  
}