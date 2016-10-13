package be.cmpg.cancer.pcawg

@Deprecated
object FunSeq2Bootstrap {

  /*
  val parser = FunSeq2Helper.getBaseArgParser("SSA. FunSeq analysis")
  
  
  val configOpt = 
  
  if (configOpt.isEmpty) {
    println("Error in parameters.")
    parser.showUsageAsError
    System.exit(1)
  }
  
  /**
   * Funtions
   */
  
  override def getArgParser  = {
    val parser = super.getArgParser;
      parser.opt[Int]("bootstraapExperiments") action { (x, c) =>
        c.copy(bootstraapExperiments = x)
      } text ("Number of calculations with random inputs to calculate the bootstraap support. Large values increase processing time. (1000 by default)")
      
      parser.opt[Double]("minBootstraapSupport") action { (x, c) =>
        c.copy(minBootstraapSupport = x)
      } text ("Minimum bootstrap support to be included in the network (0.95 by default)")
      
       parser.opt[Boolean]("useCGC") action { (x, c) =>
        c.copy(useCGC = x)
      } text ("Use the Census of Cancer Genes in the True Positives set (true by default)")
    
        parser.opt[Double]("ppv") action { (x, c) =>
        c.copy(ppv = x)
      } text ("PPV to find the closest peak.")
      parser
  }
  
  def run(config:Config, seedGenes:Set[Gene], network:Network, geneFunSeqScore:List[FunSeqData]) = {
 
 *  
 */
}