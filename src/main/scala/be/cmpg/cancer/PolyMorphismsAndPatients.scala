package be.cmpg.cancer

import be.cmpg.graph.Gene

class PolyMorphismsAndPatients {
  
}

case class Annotation (val geneSymbol: String, val impact: String)
case class Polymorphism (val geneSymbol: String, val source: String = "any", val score: Double = 1, _type: String = "any")
case class PolymorphismKey (val gene: Gene, val sample: String, val location: Double)
case class SampleInfo (val sample: String, val subtype: String)

case class GeneAnnotation(gene: Gene, chrm: String, start: Int, end: Int) {
  lazy val size = end - start
}
case class AnnotationInfo(closeGenes: Seq[GeneAnnotation], sumGenesSize: Int)