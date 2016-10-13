package be.cmpg.cancer.pcawg

@Deprecated
object FunSeq2Helper {

  def getBaseArgParser(name:String) = new scopt.OptionParser[Config](name) {
    head("Mutual exclusivity analysis by Small Subnetwork Analysis", "0.1")
    opt[Int]('i', "iterations") action { (x, c) =>
      c.copy(iterations = x)
    } text ("iterations is the number of iterations (default: 5000)")

    opt[Double]('r', "reinforcement") action { (x, c) =>
      c.copy(reinforcement = x)
    } text ("The reinforncement for succesful learning (default: 0.0005)")

    opt[Double]('f', "forgetfulness") action { (x, c) =>
      c.copy(forgetfulness = x)
    } text ("The forgetfulness after each iteration (default: 0.9995)")
    
    opt[Double]('c', "convergence") action { (x, c) =>
      c.copy(convergence = x)
    } text ("The convergence required calculated as the RMSE between the posterior distributions (default: 0.0, i.e. finish by maxing iterations)")

    opt[String]('o', "outputPrefix") action { (x, c) =>
      c.copy(outputPrefix = x)
    } text ("The prefix for the analisis output files (default: ME)")

    opt[Seq[String]]('n', "refNetwork") action { (x, c) =>
      c.copy(refNetwork = x)
    } text ("The reference network (biogrid, kinase-substrate, encode or/and files. biogrid, kinase-substrate and encode used by default.)")

    opt[Boolean]('u', "useRank") action { (x, c) =>
      c.copy(useRank = x)
    } text ("Rank the gene scores instead of the scaled values (default: true)")

    opt[String]("dir") required() action { (x, c) =>
      c.copy(dir = x)
    } text ("File Directory")
    
    opt[String]("samples") action { (x, c) =>
      c.copy(samplesFile = Some(x))
    } text ("File with the samples to analyse (all samples by default)")
    
    /*
    opt[Boolean]("codingMutations") required() action { (x, c) =>
      c.copy(codingMutations = x)
    } text ("True: Coding Mutations. False: Non Coding Mutations")
     */
    
    
    opt[Seq[String]]("mutType") action { (x, c) =>
      c.copy(mutType = x)
    } text ("annotations: Intron,Promoter,UTR,missense_variant,splice_variant,stop_gained,stop_lost (All by default)")
    
    opt[Int]("outputGenes") action { (x, c) =>
      c.copy(outputGenes = x)
    } text ("Number of genes to select (default: 200)")
    
    opt[Double]("funSeq2Threshold") action { (x, c) =>
      c.copy(funSeq2Threshold = x)
    } text ("Threhold to filter FunSeq2 scores (default: 2.0)")
    
    
    opt[Boolean]("printInputs") action { (x, c) =>
      c.copy(printInputs = x)
    } text ("Print a file with Sample Gene Score with the parsed input (only works when one mutation type is used as input)")
    
    opt[Seq[String]]('d', "debug") action { (x, c) =>
      c.copy(debug = x)
    } text ("Debug mode. List of genes that want to be observed (e.g. -d TP53,MYC,")

    
    help("help")
    
  }
  
}

case class Config (
    iterations: Int = 5000,
    reinforcement: Double = 0.0005,
    forgetfulness: Double = 0.9995,
    refNetwork: Seq[String] = Seq("biogrid", "kinase-substrate", "encode"),
    useRank: Boolean = true,
    outputPrefix: String = "PAN_FunSeq2",
    convergence: Double = 0.0,
    outputGenes:Int = 200,
    dir:String = "mumbai/",
    //codingMutations: Boolean = true, // misssense, stop_codon, etc. vs NCD (coding vs non-coding)
    funSeq2Threshold:Double = 2.0,
    mutType: Seq[String] = List("Intron", "Promoter", "UTR", "missense_variant", "splice_variant", "stop_gained", "stop_lost"), // Only for NCD: enhancer, promoter, lncrna 
    printInputs:Boolean = false,
    debug:Seq[String]=List(),
    samplesFile:Option[String]=None,
    
    /*
     * Bootstraap params
     */
    bootstraapExperiments:Int = 1000,
    minBootstraapSupport:Double = 0.95,
    useCGC:Boolean=true,
    ppv:Double=0.6
    )