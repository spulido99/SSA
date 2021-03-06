package be.cmpg.cancer

import scopt.OptionParser
import java.io.File
import scala.collection.mutable.HashMap

object MutualExclusivityPrepareInput extends App {

  case class InputConfig(
    maf: Option[File] = None,
    expression: Map[String, File] = Map(),
    seedGenesMutations: Int = 1,
    acceptedCorrelationQVal: Double = 0.05,
    acceptedGisticThreshold: Int = 2,
    maxQtyMutations: Int = 500,
    maxQtyCopyNumber: Int = 500,
    output: String = "SSAME_input")

  val helper = new CancerHelper
    
  val parser = new scopt.OptionParser[InputConfig]("SSA.ME.") {
    opt[String]('o', "outputPrefix") action { (x, c) =>
      c.copy(output = x)
    } text ("The name to be used in the output files (XXX.m2 and XXX.glst).")
    
    opt[Int]('s', "seedGenesMutations") action { (x, c) =>
      c.copy(seedGenesMutations = x)
    } text ("The number of mutated samples required in a gene to be included (default: 1)")
    
    opt[Int]("maxQtyMutations") action { (x, c) =>
      c.copy(maxQtyMutations = x)
    } text ("Maximum number of mutations for a sample. If more, the sample is removed. (default: 500)")
    
    opt[Int]("maxQtyCopyNumber") action { (x, c) =>
      c.copy(maxQtyCopyNumber = x)
    } text ("Maximum number of copy number variations for a sample. If more, the sample is removed. (default: 500)")
    
    opt[File]('m', "maf") action { (x, c) =>
      c.copy(maf = Some(x))
    } text ("Mutation .maf file")
    opt[Map[String, File]]('e', "expression") action { (x, c) =>
      c.copy(expression = x)
    } text ("expression file and GISTIC files (... -e corr=<file1>,gistic=<file2> OR" +
        " -e gistic=<file1>,exp=<file2>,cnv_thresholds=<file3> ...)." +
        " gistic folder should be a GISTIC output folder containing the files: all_thresholded.by_genes.txt, amp_genes.conf_99.txt and del_genes.conf_99.txt. " +
        " Correlation file should be tab delimited file (gene, corr, p-value, q-value).")
    opt[Double]("acceptedCorrelationQVal") action { (x, c) =>
        c.copy(acceptedCorrelationQVal = x)
    } text ("q-value of the correlation to take into account (default: 0.05)")
    
    opt[Double]("acceptedGisticThreshold") action { (x, c) =>
        c.copy(acceptedCorrelationQVal = x)
    } text ("Absolute value threshold to accept a gistic call using all_thresholded.by_genes.txt file (default: 2)")
  }

  parser.parse(args, InputConfig()) match {

    case Some(config) =>
       val genePatientMatrix = {

        //val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]  
        val mutationMatrix = if (config.maf.isDefined) {
        	println("Loading mutation file...")        
        	helper.loadMaf(config.maf.get, maxQtyMutations = config.maxQtyMutations)
        } else {
          Map[PolimorphismKey, Polimorphism]()
        }
        
        val copyNumberMatrix = if (config.expression.size >= 2) {
          println("Loading expression...")
          helper.loadExpression(config.expression, config.acceptedCorrelationQVal, config.acceptedGisticThreshold, config.maxQtyCopyNumber)
        } else {
          Map[PolimorphismKey, Polimorphism]()
        }
        
        copyNumberMatrix ++ mutationMatrix 
        
        //genePatientMatrix.toMap
      }
       
       helper.printMutationMatrixFiles(config.output, config.seedGenesMutations, genePatientMatrix)
       println("Mutation matrix file printed: "+config.output)        
      
    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }
}
