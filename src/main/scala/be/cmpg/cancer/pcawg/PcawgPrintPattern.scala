package be.cmpg.cancer.pcawg

import be.cmpg.cancer.CancerHelper
import scala.io.Source
import be.cmpg.graph.Gene
import be.cmpg.cancer.MutualExclusivityPrintPattern
import org.json.JSONArray
import be.cmpg.cancer.ColorInfo
import be.cmpg.cancer.MutualExclusivityNetworkManager

object PcawgPrintPattern extends App {

  val helper = new CancerHelper

  val parser = helper.getBasicArgParser("SSA.ME Print patterns from a list of genes")
  
  
  val mutType: Seq[String] = List("Intron", "Promoter", "UTR", "missense_variant", "splice_variant", "stop_gained", "stop_lost") // Only for NCD: enhancer, promoter, lncrna 
    
  
  parser.opt[String]("geneList") action { (x, c) =>
    c.copy(genes = x.split(",").toList)
  } text ("Gene list to print network and pattern (comma separated).")
  
  val configOpt = parser.parse(args, be.cmpg.cancer.Config(

    /*
     * Default Values
     */
    iterations = 1000,
    reinforcement = 0.005,
    forgetfulness = 0.996))
 
  if (configOpt.isEmpty) {
    parser.showUsageAsError
    System.exit(1)
  }
  
  val config = configOpt.get
  
  val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config)

  val geneFunSeqScore = Source.fromFile(config.input.get.replaceAll(".m2", ".tbs"))
                            .getLines()
                            .map { line => line.split("\t") }
                            .map { f => FunSeqData(f(0), Gene(f(1)), 0.0, f(2))}
                            .toList
  
  
  val rankedGenes = config.genes.map { Gene(_) }
                            
  val byGene = geneFunSeqScore
                  .filter { d => rankedGenes.contains(d.gene) }
                  .groupBy { g => g.gene }
                  
                  
  val exitName = config.outputPrefix+".pcawg"
                  
  val pie: Map[Gene, Map[String, Any]] = rankedGenes.map { g =>
                     val l = byGene.getOrElse(g, List())
                     val byType = l.groupBy { _.typ }.mapValues { _.size }
                     val info = mutType.map { typ => byType.getOrElse(typ, 0) }
                     (g, Map("mutTypes" -> new JSONArray(info.toArray)))
                  }.toMap
                  
  
  val networkManager = new MutualExclusivityNetworkManager(
    network = network,
    genePatientMatrix = genePatientMatrix,
    minimumSamplesAltered = 0,
    pheromone = 0,
    evaporation = 0,
    ranked = false)
  
  MutualExclusivityPrintPattern.printPattern(exitName, rankedGenes, networkManager, Map(), pie, colorInfo=ColorInfo("mutTypes", new JSONArray(mutType.toArray)))
  
  println("Network in " + exitName + "_network.html ["+rankedGenes.size+" genes selected]")
  

  
}