package be.cmpg.cancer

import be.cmpg.graph.Interaction
import be.cmpg.walk.SubNetworkSelector
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import au.com.bytecode.opencsv.CSVReader
import java.util.concurrent.Callable
import scala.collection.JavaConversions.asScalaBuffer
import scala.collection.JavaConversions.mutableSeqAsJavaList
import java.io.FileReader
import be.cmpg.graph.Gene

object CalculateSubnetworkScore extends App {

  val genes = List("PTEN","PIK3CA","CDH2","CDH1","VAV2","TTN","NCKAP5","TP53","PLK3","PIK3R1")
  
  val genePatientMatrix = {

    println("Loading mutation matrix file...")
    val genePatientMatrix =
      this.loadMutationMatrixFiles("SSAME_input.m2")

    genePatientMatrix.toMap
  }

  val mutationsPerGene = genePatientMatrix.groupBy(_._1.gene.name).map(e => (e._1, e._2.map(e => (e._1.sample, e._2))))

  val mutatedSamplesByGene = mutationsPerGene.keys.map(g => (g, mutationsPerGene(g).keySet)).toMap

  def loadMutationMatrixFiles(input: String): HashMap[PolimorphismKey, Polimorphism] = {

    val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]
    new CSVReader(new FileReader(input), '\t').readAll().foreach { fields =>
      val sample = fields(0)
      fields.drop(1).foreach { gene =>
        // Location not yet implemented here !
        val key = PolimorphismKey(Gene(gene), sample, 0)
        val value = Polimorphism(gene)

        genePatientMatrix.put(key, value)
      }
    }

    genePatientMatrix
  }
  
 println (scoreSubnetwork(genes))

  def scoreSubnetwork(subnetwork: List[String], selector: Option[SubNetworkSelector] = None): Double = {
    val allGenes = subnetwork.toSet

    val genes = allGenes.filter { g =>
      {
        val mutationsInGene = mutationsPerGene.getOrElse(g, Map()).size
        (1 < mutationsInGene && mutationsInGene < 500)
      }
    }

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
    val orderedGenes = genes.toList.sortBy(g => -1 * mutatedSamplesByGene(g).size)

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

    /*
       * Similar to Dendrix
       */
    // Make a list of samples to be analysed (samples with mutations in those genes)
    val pendingSamples = new HashSet ++ genes.map { g => mutatedSamplesByGene(g) }.reduceLeft(_ ++ _)

    val scorePerGene = orderedGenes.map {
      gene =>
        {
          var geneScore = 0.0
          val usedSamples = new HashSet[String]
          for (sample <- pendingSamples) {
            // Mutated samples should be 1 for perfect mutual exclusivity
            // The number samples that have the gene mutated
            if (mutatedSamplesByGene(gene).contains(sample)) {
              val mutatedSamples = orderedGenes.map { otherGene => if (mutatedSamplesByGene(otherGene).contains(sample)) 1.0 else 0.0 }.reduceLeft(_ + _)
              usedSamples += sample
              //geneScore += (genes.size - mutatedSamples) / genes.size
              geneScore += 1.0 / mutatedSamples
            }
          }
          pendingSamples --= usedSamples
          (gene, geneScore)
        }
    }
      .toMap
    return scorePerGene.map { sg => math.sqrt(sg._2) }.sum / allGenes.size
  }

}