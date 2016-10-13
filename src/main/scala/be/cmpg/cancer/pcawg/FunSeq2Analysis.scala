package be.cmpg.cancer.pcawg

import scala.io.Source
import be.cmpg.cancer.CancerHelper
import be.cmpg.graph.Network
import be.cmpg.graph.Gene
import be.cmpg.cancer.MutualExclusivityPrintPattern
import scala.util.Try
import scala.util.Failure
import scala.util.Success
import java.io.PrintWriter
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import be.cmpg.graph.Gene
import org.json.JSONArray
import be.cmpg.cancer.ColorInfo

@Deprecated
object FunSeq2Analysis extends App {

  val parser = getArgParser
  
  val helper = new CancerHelper

  val configOpt = parser.parse(args, Config())
  
  if (!configOpt.isDefined) {
    parser.showUsageAsError
    System.exit(1)
  }
  
  val config = configOpt.get

  val badGenes = Set("UBC", "TTN").map(Gene(_))
  val interactions = helper.loadNetwork(config.refNetwork).filter { i => badGenes.intersect(i.genes).size == 0 }
  val network = new Network(interactions)
  
  println("Interactions: "+interactions.size)
  
//  var samples = if (config.samplesFile.isDefined) Source.fromFile(config.samplesFile.get, "latin1") else None
  val samples = if (config.samplesFile.isDefined) Source.fromFile(config.samplesFile.get, "latin1").getLines.map { line => line.split(",") }.flatten.toSet else Set[String]()
  
  /*
   * Read all mutations files: From chr1 to chr23 and X, Y
   * 
   * - Try to Read tab separated line
   * - Filter short, empty or incomplete lines
   * - map to data structure
   * - filter by coding or noncoding mutations
   * - get gene name and fnuseq score depending if coding or non coding
   * - filter by genes in the network
   */
  
  
  val geneFunSeqScore = Source.fromFile(config.dir + "mutations.funseq2.tbs")
                            .getLines()
                            .map { line => line.split("\t") }
                            .map { f => FunSeqData(f(0), Gene(f(1)), f(2).toDouble, f(3))}
                            .filter { d => config.mutType.contains(d.typ) }
                            .filter { d => d.score > config.funSeq2Threshold }
                            .toList
    
    /*(List.range(1, 23) ++ List("X", "Y")).toStream.map { input => 
                            val filename = config.dir+"chr"+input+".vcf.data.tsv"
                            Try {
                              Source.fromFile(filename, "latin1")
                              .getLines
                              .map { line => line.split("\t") }
                              .filter { _.size >= 7 }
                              .map { f => Data(f(0), f(1), f(2).toInt, f(3), f(4), if (f(6).isEmpty()) 0.0 else f(6).toDouble, if(f.size == 8 && !f(7).isEmpty()) f(7).split(":")(0).toDouble else 0.0)}                            
                              .filter { d => 
                                  val coding = config.codingMutations && config.mutType.get.contains(d.muttyp)
                                  val nonCoding = !config.codingMutations && !config.mutType.get.filter { mutType => d.gene.indexOf(mutType) > 0 }.isEmpty
                                  if (config.codingMutations) coding else nonCoding
                              }.map { d => 
                                val gene  = if (config.codingMutations) d.gene     
                                            else {
                                                	d.gene.split(",").find { geneAnn => 
                                                    !config.mutType.get.filter { mutType => geneAnn.indexOf(mutType) > 0 }.isEmpty 
                                                  }.getOrElse("").split("\\(")(0)
                                              }
                                val score = if (config.codingMutations) d.cdsScore else d.ncdScore
                                FunSeqData(Gene(gene), d.sample, score)
                              }.filter { fsd => network.genes contains fsd.gene }
                              .filter { fsd => (samples.isEmpty || samples.contains(fsd.sample)) } 
                            } match {
                              case Success(v) =>
                                println("File loaded: "+filename)
                                Some(v)
                              case Failure(e) => 
                                println("*ERROR* File not loaded: "+filename)
                                None
                            }
                        }.flatten.flatten.toList
                        * 
                        
  
  if (config.printInputs) {
    
    val writer = new PrintWriter("ssame.funseq2.tsv", "UTF-8");
    geneFunSeqScore.foreach { fsd =>
      writer.println(fsd.sample+"\t"+fsd.gene.name+"\t"+fsd.score)
    }
    writer.close
    System.exit(0);
  }
  */
  
  /*
   * Some statas about the FunSeqScores
   */
  val scoreThreshold = {
    val cgcStats = new DescriptiveStatistics
    val ncgStats = new DescriptiveStatistics
    val allStats = new DescriptiveStatistics
    
    geneFunSeqScore.groupBy { _.gene }.foreach { case (g, values) =>
      val scoresum = values.foldLeft(0.0)(_+_.score) 
      
      if (helper.cgc contains g) {
        cgcStats.addValue(scoresum)
      }
      if (helper.ncg contains g) {
    	  ncgStats.addValue(scoresum)
      }
      
      allStats.addValue(scoresum)
    }
    
    println("CGC genes FunSeq2 score stats")
    println(cgcStats)
    println("NCG genes FunSeq2 score stats")
    println(ncgStats)
    println("All genes FunSeq2 score stats")
    println(allStats)

    math.min(cgcStats.getPercentile(50), ncgStats.getPercentile(50))
  };
  println("Threshold to use: "+scoreThreshold)
                        
  
  /*
   * Seed genes will be genes with a sum of FunSeq score > 50
   */
  val seedGenes = geneFunSeqScore.groupBy { _.gene }.mapValues(_.foldLeft(0.0)(_+_.score)).filter(_._2 > scoreThreshold).keySet
  println("Seed genes: "+seedGenes.size)
  
  val (rankedGenes, networkManager) = run(config, seedGenes, network, geneFunSeqScore)
  
  val exitName = config.outputPrefix+".funseq"
  
  val byGene = geneFunSeqScore
                  .filter { d => rankedGenes.contains(d.gene) }
                  .groupBy { g => g.gene }
                  
  val pie: Map[Gene, Map[String, Any]] = rankedGenes.map { g =>
                     val l = byGene.getOrElse(g, List())
                     val byType = l.groupBy { _.typ }.mapValues { _.size }
                     val info = config.mutType.map { typ => byType.getOrElse(typ, 0) }
                     (g, Map("mutTypes" -> new JSONArray(info.toArray)))
                  }.toMap
  
  MutualExclusivityPrintPattern.printPattern(exitName, rankedGenes, networkManager, Map(), pie, colorInfo=ColorInfo("mutTypes", new JSONArray(config.mutType.toArray)))
  
  println("Network in " + exitName + "_network.html ["+rankedGenes.size+" genes selected]")
  
  
  /**
   * Funtions
   */
  
  def getArgParser = FunSeq2Helper.getBaseArgParser("SSA. FunSeq analysis")
  
  def run(config:Config, seedGenes:Set[Gene], network:Network, geneFunSeqScore:List[FunSeqData]) = {
    val networkManager = new FunSeq2NetworkManager(
      network = network,
      geneFunSeqScore = geneFunSeqScore,
      pheromone = config.reinforcement,
      evaporation = config.forgetfulness,
      ranked = config.useRank,
      convergenceThreshold = config.convergence)
  
    if (!config.debug.isEmpty)
          networkManager.debug = Some(config.debug.map { Gene(_) }.toSet, 20)
                          
    networkManager.run(config.iterations, seedGenes, Runtime.getRuntime().availableProcessors(), storeNodePHistory = true)
    
    val rankedGenes = if (config.outputGenes == 0) networkManager.rankedGenes.toList else networkManager.getRankedAllGenes().take(config.outputGenes).toList
 
    (rankedGenes, networkManager)
  }
  
}

case class FunSeqData(sample:String, gene:Gene, score:Double, typ:String)
case class Data(sample:String, chromosome:String, location:Int, gene:String, muttyp:String, cdsScore:Double, ncdScore:Double)
