package be.cmpg.cluster

import scopt.OptionParser

object ArgumentsParser {

  def getBasicArgParser(name: String): OptionParser[PvalueConfig] = {

    new scopt.OptionParser[PvalueConfig](name) {

      head("Small subnetwork analysis of p-values", "0.1")

      opt[Int]('i', "iterations") action { (x, c) =>
        c.copy(iterations = x)
      } text ("iterations is the number of iterations (default: 50000)")

      opt[Double]('r', "reinforcement") action { (x, c) =>
        c.copy(reinforcement = x)
      } text ("The reinforncement for succesful learning (default: 0.0005)")

      opt[Double]('f', "forgetfulness") action { (x, c) =>
        c.copy(forgetfulness = x)
      } text ("The forgetfulness after each iteration (default: 0.9996)")

      opt[String]('o', "outputFolder") required () action { (x, c) =>
        c.copy(outputFolder = x)
      } text ("The folder to where the output should be written.")

      opt[Double]('s', "seedGenesMutations") action { (x, c) =>
        c.copy(seedGenesMutations = x)
      } text ("The minimum -log(p) value for a gene to be considered a seed. Default is 0.1")

      opt[String]('p', "outputPrefix") action { (x, c) =>
        c.copy(outputPrefix = x)
      } text ("The prefix for the analisis output files (default: CLUST)")

      opt[String]('n', "refNetwork") required () action { (x, c) =>
        c.copy(refNetwork = x)
      } text ("Path to the network file. Interactions must be undirected. Directed interactions are NOT supported yet.")

      opt[Boolean]("patterns") action { (x, c) =>
        c.copy(patterns = x)
      } text ("calculate 5 best subnetworks of every gene as ME pattern for that gene. Default is false")

      opt[Boolean]('u', "useRank") action { (x, c) =>
        c.copy(useRank = x)
      } text ("Rank the gene scores instead of the scaled values (default: false)")

      opt[String]('m', "input_list") required () action { (x, c) =>
        c.copy(input_list = x)
      } text ("Path to the input list containing genes with their p-values in CSV format.")

      opt[Int]('p', "outputGenes") action { (x, c) =>
        c.copy(outputGenes = x)
      } text ("Number of nodes to be published to the html network output (default: 100).")

      opt[Int]("processors") action { (x, c) =>
        c.copy(processors = x)
      } text ("Number parallel processors to use (default: half of available processors).")

      opt[String]('a', "inputFolder") required () action { (x, c) =>
        c.copy(inputFolder = x)
      } text ("Path to folder where network html scaffold is located. SHOULD BE REFACTORED.")

      opt[String]('e', "excludedGenes") action { (x, c) =>
        c.copy(excludedGenes = x)
      } text ("Path to file which is a list with genes to exclude from the analysis (e.g. large hubs known to influence the analysis too much).")
    }
  }
}
case class PvalueConfig(

  iterations: Int = 5000,
  reinforcement: Double = 0.0005,
  forgetfulness: Double = 0.9996,
  outputFolder: String = "",
  seedGenesMutations: Double = 0.1,
  outputPrefix: String = "CLUST",
  refNetwork: String = "",
  patterns: Boolean = false,
  useRank: Boolean = false,
  input_list: String = "",
  outputGenes: Int = 100,
  genes: List[String] = List(),
  inputFolder: String = "",
  excludedGenes: String = "",
  processors: Int = Runtime.getRuntime().availableProcessors() / 2)
