package be.cmpg.cancer.ME

import scopt.OptionParser
import java.io.File
import scala.collection.mutable.HashMap
import scala.io.Source
import be.cmpg.cancer.Polymorphism
import be.cmpg.cancer.PolymorphismKey
import be.cmpg.cancer.DataLoader
import be.cmpg.cancer.PrintFiles

object MutualExclusivityPrepareInput extends App {

  case class InputConfig(
    maf: Option[File] = None,
    expression: Map[String, File] = Map(),
    seedGenesMutations: Int = 1,
    acceptedCorrelationQVal: Double = 0.05,
    acceptedGisticThreshold: Int = 2,
    maxQtyMutations: Int = 500,
    maxQtyCopyNumber: Int = 500,
    GenesToBeConsidered: String = "",
    outputFolder:String = "",
    inputFolder:String="",
    output: String = "SSAME_input")

  val parser = new scopt.OptionParser[InputConfig]("SSA.ME.") {
    opt[String]('p', "outputPrefix") action { (x, c) =>
      c.copy(output = x)
    } text ("The name to be used in the output files (XXX.m2 and XXX.glst).")

    opt[String]('o', "outputFolder") required() action { (x, c) =>
      c.copy(outputFolder = x)
    } text ("The folder to write the output to.")

    opt[Int]('s', "seedGenesMutations") action { (x, c) =>
      c.copy(seedGenesMutations = x)
    } text ("The number of mutated samples required in a gene to be included (default: 1)")

    opt[String]('l', "GenesToBeConsidered") action { (x, c) =>
      c.copy(GenesToBeConsidered = x)
    } text ("A list of genes which genes should only be considered. If left empty all genes are considered")

    opt[Int]("maxQtyMutations") action { (x, c) =>
      c.copy(maxQtyMutations = x)
    } text ("Maximum number of mutations for a sample. If more, the sample is removed. (default: 500)")

    opt[String] ('i',"inputFolder") required() action { (x, c) =>
      c.copy(inputFolder = x)
    } text ("Folder where the input files are situated (network files, possibly known cancer gene files)")
    
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

    opt[Int]("acceptedGisticThreshold") action { (x, c) =>
      c.copy(acceptedGisticThreshold = x)
    } text ("Absolute value threshold to accept a gistic call using all_thresholded.by_genes.txt file (default: 2)")
  }

  parser.parse(args, InputConfig()) match {

    case Some(config) =>
      val genePatientMatrix = {

        val GenesToBeConsidered = if (config.GenesToBeConsidered == "") { List() } else { Source.fromFile(config.GenesToBeConsidered).getLines.drop(1).map(line => line.split("\t")(0)).toList }
        //val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]  
        val mutationMatrix = if (config.maf.isDefined) {
          println("Loading mutation file...")
          DataLoader.loadMaf(config.maf.get, maxQtyMutations = config.maxQtyMutations, GenesConsidered = GenesToBeConsidered)
        } else {
          Map[PolymorphismKey, Polymorphism]()
        }

        val copyNumberMatrix = if (config.expression.size >= 2) {
          println("Loading expression...")
          DataLoader.loadExpression(config.expression, config.inputFolder, config.acceptedCorrelationQVal, config.acceptedGisticThreshold, config.maxQtyCopyNumber, GenesConsidered = GenesToBeConsidered)
        } else {
          Map[PolymorphismKey, Polymorphism]()
        }

        copyNumberMatrix ++ mutationMatrix

        //genePatientMatrix.toMap
      }

      PrintFiles.printMutationMatrixFiles(config.output, config.seedGenesMutations, genePatientMatrix, false,config.outputFolder)
      println("Mutation matrix file printed: " + config.output)

    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }
}