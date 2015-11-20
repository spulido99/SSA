package be.cmpg.qtl

/**
 * The PPI network should be completed with marker information.
 * 
 * Edges between each gene and the marker sites that happen in that 
 * gene or its regulatory region. The marker sites are joined by edges
 * based on its spatial chromosomal position.
 */
object QtlEpistasisAnalysis extends App {

  // Marker NODES (Gene) should have data with "chromosome" and "pos"
  
}

object QtlConstants {
  val chromosome = "chr"
  val position = "pos"
    
  val startLoci = "startLoci"
  val endLoci = "endLoci"
}