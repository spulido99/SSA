package be.cmpg.cancer
import scopt.OptionParser
import java.io.File
import scala.io.Source

// This class contains arguments and default values needed for all cancer projects of SSA.
object ArgumentsParser {
  def getBasicArgParser(name: String): OptionParser[Config] = {

    new scopt.OptionParser[Config](name) {

      head("Mutual exclusivity analysis by Small Subnetwork Analysis", "0.1")
      opt[Int]('i', "iterations") action { (x, c) =>
        c.copy(iterations = x)
      } text ("iterations is the number of iterations (default: 50000)")

      opt[Double]('r', "reinforcement") action { (x, c) =>
        c.copy(reinforcement = x)
      } text ("The reinforncement for succesful learning (default: 0.0005)")

      opt[String]('e', "excluded-genes") action {(x, c) =>
        c.copy(excludedGenes = x)
      } text ("Path to a file containing the excluded genes (one gene name per line).")
      
      opt[String]('o', "outputFolder") action { (x, c) =>
        c.copy(outputFolder = x)
      } text ("The folder to where the output should be written.")

      opt[Double]('f', "forgetfulness") action { (x, c) =>
        c.copy(forgetfulness = x)
      } text ("The forgetfulness after each iteration (default: 0.9996)")

      opt[Int]('s', "seedGenesMutations") action { (x, c) =>
        c.copy(seedGenesMutations = x)
      } text ("The number of mutated samples required in a gene to be a seed (default: 1)")

      opt[String]('p', "outputPrefix") action { (x, c) =>
        c.copy(outputPrefix = x)
      } text ("The prefix for the analisis output files (default: ME)")

      opt[Seq[String]]('n', "refNetwork") action { (x, c) =>
        c.copy(refNetwork = x)
      } text ("The reference network as a sequence of strings ('HT' (hint), 'hiII14' (interactome), 'reactome' or 'MEMo' (MEMo reference network)). HT, hiII14 and reactome together used by default. note that interactions are viewed as undirected interactions.)")

      opt[Boolean]("patterns") action { (x, c) =>
        c.copy(patterns = x)
      } text ("calculate 5 best subnetworks of every gene as ME pattern for that gene. Default is false")

      opt[Boolean]('u', "useRank") action { (x, c) =>
        c.copy(useRank = x)
      } text ("Rank the gene scores instead of the scaled values (default: true)")

      opt[Boolean]('h', "statistical") action { (x, c) =>
        c.copy(statistical = x)
      } text ("Use an hypergeometric test to search for mutual exclusivity. Will not find low mutations genes. (default: false)")

      opt[Int]('g', "minMutPerGene") action { (x, c) =>
        c.copy(minMutPerGene = x)
      } text ("Ignore genes mutated in less than 'g' samples (default: 1)")

      opt[String]('a', "otherGeneList") action { (x, c) =>
        c.copy(otherGeneList = x)
      } text ("Path to file containing genes which will be included in the visualization")

      opt[String]('m', "input_matrix") required () action { (x, c) =>
        var inputMatrix = x
        if (!x.endsWith(".m2"))
          inputMatrix = x + ".m2"
        c.copy(input = Some(inputMatrix))
      } text ("Mutation matrix .m2 file (created with option ME_input )")

      opt[String]('b', "inputFolder") required () action { (x, c) =>
        c.copy(inputFolder = x)
      } text ("Path to the folder containing the networks and the cgc/ncg/other genes.")

      opt[Int]('p', "outputGenes") action { (x, c) =>
        c.copy(outputGenes = x)
      } text ("Number of nodes to be published to the html network output (default: 500).")

      opt[Int]("processors") action { (x, c) =>
        c.copy(processors = x)
      } text ("Number parallel processors to use (default: half of available processors).")

      opt[File]("subtypeFile") action { (x, c) =>
        c.copy(subtypeFile = Some(x))
      } text ("Samples per subtype (Tab Delimited File with 2 columns: Sample Subtype")

      opt[String]("subtype") action { (x, c) =>
        c.copy(subtype = Some(x))
      } text ("subtype to analyse")

      opt[Seq[String]]('d', "debug") action { (x, c) =>
        c.copy(debug = Some(x))
      } text ("Debug mode. List of genes that want to be observed (e.g. -d TP53,MYC,")
      
      help("help")
    }
  }
}
case class Config(
  /*
   * Base parameters
   */
  iterations: Int,
  reinforcement: Double,
  forgetfulness: Double,
  refNetwork: Seq[String] = Seq("HT", "hiII14", "reactome"),
  useRank: Boolean = true,
  statistical: Boolean = false,
  minMutPerGene: Int = 1,
  inputFolder: String = "",
  seedGenesMutations: Int = 1,
  outputPrefix: String = "ME",
  input: Option[String] = None,
  outputGenes: Int = 500,
  otherGeneList: String = "",
  subtype: Option[String] = None,
  subtypeFile: Option[File] = None,
  debug: Option[Seq[String]] = None,
  outputFolder: String = "",
  excludedGenes:String ="",
  patterns: Boolean = false,
  processors: Int = Runtime.getRuntime().availableProcessors() / 2,

  /*
   * P-value calculation parameters
   */
  pvalueExperiments: Int = 1000,
  closestGenes: Int = 5,
  genomicdistnace: Int = 50000,
  pvaluetype: String = "distance",

  /*
   * Bootstrap calculation parameters
   */
  bootstrapExperiments: Int = 100,
  minBootstrapSupport: Double = 0.9,
  useCGC: Boolean = true,
  useNCG: Boolean = true,
  bootstrapOutPutGenes: Int = 100,

  /*
   * Print patern additional parameters
   */
  genes: List[String] = List())
