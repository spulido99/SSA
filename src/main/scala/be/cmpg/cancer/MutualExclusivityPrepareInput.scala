package be.cmpg.cancer

import scopt.OptionParser
import java.io.File
import scala.collection.mutable.HashMap

object MutualExclusivityPrepareInput extends App {

  case class InputConfig(
    maf: Option[File] = None,
    expression: Map[String, File] = Map(),
    seedGenesMutations: Int = 1,
    output: String = "SSAME_input")

  val helper = new CancerHelper
    
  val parser = new scopt.OptionParser[InputConfig]("SSA.ME.") {
    opt[String]('o', "outputPrefix") action { (x, c) =>
      c.copy(output = x)
    } text ("The name to be used in the output files (XXX.m2 and XXX.glst).")
    opt[Int]('s', "seedGenesMutations") action { (x, c) =>
      c.copy(seedGenesMutations = x)
    } text ("The number of mutated samples required in a gene to be included (default: 1)")
    opt[File]('m', "maf") required () action { (x, c) =>
      c.copy(maf = Some(x))
    } text ("Mutation .maf file")
    opt[Map[String, File]]('e', "expression") action { (x, c) =>
      c.copy(expression = x)
    } text ("expression file and GISTIC files (... -e cnv_peaks=<file1>,exp=<file2>,cnv_thresholds=<file3> ...).")
  }

  parser.parse(args, InputConfig()) match {

    case Some(config) =>
       val genePatientMatrix = {

        val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]
        println("Loading mutation file...")        
        helper.loadMaf(genePatientMatrix, config.maf.get)
        
        if (config.expression.size == 3) {
          println("Loading expression...")
          helper.loadExpression(genePatientMatrix,
            config.expression("cnv_peaks"),
            config.expression("exp"),
            config.expression("cnv_thresholds"))
        }
        
        genePatientMatrix.toMap
      }
       
       helper.printMutationMatrixFiles(config.output, config.seedGenesMutations, genePatientMatrix)
       println("Mutation matrix file printed: "+config.output)        
      
    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }
}