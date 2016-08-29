package be.cmpg.statistical

import be.cmpg.cancer.PolymorphismKey
import be.cmpg.cancer.Polymorphism
import scala.util.Random
import java.io.File
import be.cmpg.graph.Gene
import be.cmpg.cancer.AnnotationInfo
import scala.collection.mutable.HashMap
import scala.collection.mutable.ArrayBuffer

// This object contains multiple ways to randomize gene patient matrices
object DataRandomizations {
  
    // Set up random generator
  val random = new Random(System.nanoTime())
  def getCloseGenesRandomizedGenePatientMatrix(genePatientMatrix: Map[PolymorphismKey, Polymorphism], annotationMap: Map[Gene, AnnotationInfo]): Map[PolymorphismKey, Polymorphism] = {

    val toReturn = new HashMap[PolymorphismKey, Polymorphism]

    genePatientMatrix.keys.foreach { key =>

      val annotation = annotationMap.get(key.gene)
      if (annotation.isDefined) {
        val annotationInfo = annotation.get
        var newGene: Gene = annotationInfo.closeGenes(0).gene
        var ran = random.nextInt(annotationInfo.sumGenesSize)

        annotationInfo.closeGenes
          .takeWhile(_ => ran > 0)
          .foreach { annotation =>
            ran -= annotation.size
            newGene = annotation.gene
          }
        // Location not yet implemented here !
        toReturn.put(PolymorphismKey(newGene, key.sample, 0), Polymorphism(newGene.name))
      }
    }

    toReturn.view.toMap
  }

  // Generates bootstrap samples
  def bootstrapSamples(alterationsBySample: Map[String, List[(PolymorphismKey, Polymorphism)]]): Map[PolymorphismKey, Polymorphism] = {

    val samples = alterationsBySample.keys.toArray

    List.tabulate(samples.size)(n =>

      // Location not yet implemented here !
      alterationsBySample(samples(random.nextInt(samples.size))).map(e => (PolymorphismKey(e._1.gene, "Sample" + n, 0), Polymorphism(e._2.geneSymbol)))).flatten.toMap
  }

  def randomizeGeneNames(genePatientMatrix: Map[PolymorphismKey, Polymorphism], genesSet: Set[Gene]): Map[PolymorphismKey, Polymorphism] = {

    val originalGenes = genesSet.toList
    val randomizedGenes = random.shuffle(originalGenes)

    val old2new = originalGenes.zip(randomizedGenes).toMap

    genePatientMatrix
      .filter { e => genesSet contains e._1.gene }
      .map { e =>
        val newGene = old2new(e._1.gene)
        // Location not yet implemented here !
        (PolymorphismKey(newGene, e._1.sample, 0), Polymorphism(newGene.name))
      }.toMap
  }

  def randomizeBipartiteGenePatientMatrix(genePatientMatrix: Map[PolymorphismKey, Polymorphism], samplesSet: Set[String], genesSet: Set[Gene]): Map[PolymorphismKey, Polymorphism] = {

    print("Randomizing data... ")

    val start = System.currentTimeMillis()

    val samples = samplesSet.zipWithIndex.toMap
    val genes = genesSet.zipWithIndex.toMap

    val m = Array.fill(genes.size)(Array.fill(samples.size)(false))
    genes.par.foreach { g =>
      samples.foreach { s =>
        // Location not yet implemented here !
        if (genePatientMatrix.contains(PolymorphismKey(g._1, s._1, 0))) {
          m(g._2)(s._2) = true
        }
      }
    }
    println("Done creating m. [" + (System.currentTimeMillis() - start) + " ms]")

    val Np = 100 * genePatientMatrix.size
    val d = genePatientMatrix.size / (samples.size * genes.size)

    val N = Np * math.log((1 - d) * Np / 100) / (200 * (1 - d))

    val matrix = ArrayBuffer(genePatientMatrix.toList: _*)

    println("Done creating new Matrix. [" + (System.currentTimeMillis() - start) + " ms]")

    (0 to N.toInt).par.foreach { n =>

      val r1 = random.nextInt(matrix.size)
      val r2 = random.nextInt(matrix.size)

      val ab = matrix(r1)
      val cd = matrix(r2)

      if (r1 != r2 && ab._1.gene != cd._1.gene && ab._1.sample != cd._1.sample) {

        this.synchronized {
          val a = genes(ab._1.gene)
          val b = samples(ab._1.sample)
          val c = genes(cd._1.gene)
          val d = samples(cd._1.sample)
          if (!m(a)(d) && !m(c)(b)) {
            m(a)(d) = true
            m(c)(b) = true
            m(a)(b) = false
            m(c)(d) = false
            // Location not yet implemented here !
            matrix(r1) = (PolymorphismKey(ab._1.gene, cd._1.sample, 0), ab._2)
            matrix(r2) = (PolymorphismKey(cd._1.gene, ab._1.sample, 0), cd._2)
          }
        }
      }
    }

    println("Done randomization. [" + (System.currentTimeMillis() - start) + " ms]")
    val toReturn = matrix.toMap
    println("Done. [" + (System.currentTimeMillis() - start) + " ms]")

    toReturn
  }
  
}