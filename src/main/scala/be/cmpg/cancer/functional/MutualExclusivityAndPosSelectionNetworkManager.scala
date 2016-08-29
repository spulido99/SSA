package be.cmpg.cancer.functional

import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import scala.collection.mutable.HashSet
import be.cmpg.graph.Gene
import scala.collection.Map
import scala.collection.Set
import be.cmpg.walk.SubNetworkSelector
import org.apache.commons.math3.distribution.HypergeometricDistribution
import be.cmpg.cancer.Polymorphism
import be.cmpg.cancer.PolymorphismKey

class MutualExclusivityAndPosSelectionNetworkManager(network: Network,
  genePatientMatrix: Map[PolymorphismKey, Polymorphism],// (EntrezGeneId, PatientId) -> chrom:start-end
  CADDScores:Map[Gene,Double],
  minimumSamplesAltered: Int, // minFreqGene: minimum frequency of mutation for a gene to be considered in the analysis
  maximumSamplesAltered: Int = Int.MaxValue,
  pheromone: Double = 0.05,
  evaporation: Double = 0.996,
  ranked: Boolean = false,
  minProb: Double = 0.01,
  hypergeometricTest: Boolean = false) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, ranked: Boolean, minProb: Double) {
 
  println("MEManager")
  val all_samples = genePatientMatrix.map(_._1.sample).toSet//Set[String], // all the samples names (here to not need to analyse all keys of the genePatientMatrix)
  val mutationsPerSample = genePatientMatrix.groupBy(_._1.sample).map(e => (e._1, e._2.map(e => (e._1.gene, e._2))))
  println("Samples: "+mutationsPerSample.size)
  val mutationsPerGene = genePatientMatrix.groupBy(_._1.gene).map(e => (e._1, e._2.map(e => (e._1.sample, e._2))))
  println("Genes: "+mutationsPerGene.size)
  
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

    // Get all genes which have at least 2 mutations and which do not have more mutations then some maximum (to shield from errors)
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
     * Calculate a score:
     */
    
    val mutatedSamplesByGene = genes.map( g => (g, mutationsPerGene(g).keySet)).toMap
    
    // Order genes by to analyse those with more mutations first
    val orderedGenes = genes.toList.sortBy( g => -1*mutatedSamplesByGene(g).size )
    
    // Make a list of samples to be analysed (samples with mutations in those genes)
    val pendingSamples = new HashSet ++ genes.map { g => mutatedSamplesByGene(g) }.reduceLeft(_ ++ _)
   
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

    if (hypergeometricTest) {
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
        // the score gets normalized for subnetwork size by dividing by subnetwork size which is not mentioned in the paper
       return scorePerGene.map {sg => math.sqrt(sg._2)}.sum / allGenes.size
    }
  }
}