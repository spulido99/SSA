package be.cmpg.cancer

import java.io.File
import scala.collection.mutable.HashMap
import java.io.FileWriter
import org.json.JSONObject
import org.json.JSONArray
import be.cmpg.graph.Gene
import be.cmpg.graph.Network
import be.cmpg.graph.Network
import be.cmpg.graph.interaction.NodeCostNetworkManager
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.JavaConversions._
import java.io.InputStreamReader
import scala.io.Source

object MutualExclusivityPrintPattern extends App {
  
  val helper = new CancerHelper

  val parser = helper.getBasicArgParser("SSA.ME Print patterns from a list of genes")
  
  parser.opt[String]("geneList") action { (x, c) =>
    c.copy(genes = x.split(",").toList)
  } text ("Gene list to print network and pattern (comma separated).")
  
  parser.opt[Seq[String]]("positiveGeneSetLists") action { (x, c) =>
    c.copy(positiveGeneSetLists = x)
  } text ("Files with a list of gene containing what are considered positive (COSMIC and NCG added by default). Sould be tab delimited files with the gene name in the first column.")

  val configOpt = parser.parse(args, Config(

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
  // parser.parse returns Option[C]


  val config = configOpt.get
  
  val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config)

  val networkManager = new MutualExclusivityNetworkManager(
    network = network,
    genePatientMatrix = genePatientMatrix,
    minimumSamplesAltered = 0,
    pheromone = 0,
    evaporation = 0,
    ranked = false)
  
  val otherPositiveSet = config.positiveGeneSetLists.map { geneListFile =>
                                        Source.fromFile(geneListFile, "latin1").getLines.map { _.split("\t")(0) }.map(Gene(_))
                                      }.flatten.toSet

  printPattern(config.outputPrefix+".selected", config.genes.map { n => Gene(n) }, networkManager, genePatientMatrix, otherPositiveGeneSetLists=otherPositiveSet)

  def printPattern(outputPrefix: String, genes: List[Gene], networkManager: NodeCostNetworkManager, genePatientMatrix: Map[PolimorphismKey, Polimorphism], additional: Map[Gene, Map[String, Any]] = Map(), otherPositiveGeneSetLists:Set[Gene]=Set()) = {

    val helper = new CancerHelper
                                        
    val results = new JSONObject

    val edges = new JSONArray

    val additionalNodeInfo = if (additional.isEmpty && !genePatientMatrix.isEmpty) {
      val mutationsPerGene = genePatientMatrix.groupBy(_._1.gene)

      genes.map { g =>
        val info = mutationsPerGene.get(g)
        if (info.isDefined)
          g -> Map[String, Any]("pvalue" -> math.sqrt(info.get.size), "origin" -> "coding")
        else
          g -> Map[String, Any]("pvalue" -> 0.0, "origin" -> "unkown")
      }.toMap[Gene, Map[String, Any]]
    } else {
       additional 
    }

    val out_edges = new FileWriter(outputPrefix + "_edges.tsv")
    out_edges.write("GeneA\tGeneB\tMEScore\n")
    networkManager.getAllInteractions
      .filter { i => i.to != i.from && genes.contains(i.from) && genes.contains(i.to) }
      .foreach { i =>
        {
          val score = { val score = networkManager.scoreSubnetwork(Set(i)); if (score.isNaN) 0.0 else score }

          val edgeInfo = new JSONObject()
            .accumulate("source", i.from.name)
            .accumulate("target", i.to.name)
            .accumulate("score", score)
          edges.put(edgeInfo)
          out_edges.write(i.from.name + "\t" + i.to.name + "\t" + score + "\n")
        }
      }

    out_edges.close()

    def getKnownCancerGene(gene:Gene) = {
      val toReturn = new JSONArray();
      
      toReturn.put(helper.cgc.contains(gene))
      toReturn.put(otherPositiveGeneSetLists.contains(gene))
      toReturn.put(helper.ncg.contains(gene))

      toReturn
    }  
    
    def getKnownCancerLists(gene:Gene) = {
      var toReturn = 0;
      
      if (otherPositiveGeneSetLists.contains(gene)) toReturn += 1
      if (helper.cgc.contains(gene)) toReturn += 1
      if (helper.ncg.contains(gene)) toReturn += 1
      
      toReturn
    }  
    
    val truePositiveGeneList = helper.cgc ++ helper.ncg ++ otherPositiveGeneSetLists;
    
    val nodes = new JSONArray

    val out_nodes = new FileWriter(outputPrefix + "_nodes.tsv")
    out_nodes.write("GeneSymbol\tPosteriorP\tSelected\tConvergenceIter\tKnownCancerGene\tBestSSNScore\tBestSSN\n")
    println("GeneSymbol\tOtherPositive\tCGC\tNCG")
    genes
      .foreach(g => {
        
        println(g.name+"\t"+otherPositiveGeneSetLists.contains(g)+"\t"+helper.cgc.contains(g)+"\t"+helper.ncg.contains(g))
        
        val geneInfo = new JSONObject()
          .accumulate("name", g.name)
          .accumulate("selected", genes.contains(g));
        
        geneInfo.put("knownCancerGene", getKnownCancerGene(g))
          
        val additionalInfo = additionalNodeInfo.get(g)
        if (additionalInfo.isDefined) {
          additionalInfo.get.foreach { e =>
            geneInfo.accumulate(e._1, e._2)
          }
        }
        nodes.put(geneInfo)
        try {
        	val node = networkManager.getNetwork().getNode(g)
        	out_nodes.write(g.name + "\t" + node.posteriorProbability+ "\t" + genes.contains(g) + "\t" + node.convergenceIteration + "\t" + truePositiveGeneList.contains(g) + "\t" + node.bestSubnetwork._2+ "\t" + node.bestSubnetwork._1+"\n")
        } catch {
          case t: Throwable => out_nodes.write(g.name + "\t" + "NA" + "\t" + genes.contains(g) + "\t" + "NA" + "\t" + truePositiveGeneList.contains(g) + "\t" + "NA" + "\n")
        }
      })

    out_nodes.close()

    results.put("edges", edges)
    results.put("nodes", nodes)

    val html = Source.fromInputStream(getClass.getResourceAsStream("/baseNetworkOutput.html")).mkString

    val out_json = new FileWriter(outputPrefix + "_network.html")
    out_json.write(html.replace("$${graph}$$", results.toString))

    out_json.close()

    val all_samples = genePatientMatrix.map(_._1.sample).toSet

    val geneMutScoreCounts = genes.map(gene => (gene, {
      var mutationCounts = 0
      for (sample <- all_samples) {
        if (genePatientMatrix contains (PolimorphismKey(gene, sample))) {
          mutationCounts += 1
        }
      }
      var mutualExScore = 0.0
      for (sample <- all_samples) {
        if (genePatientMatrix contains (PolimorphismKey(gene, sample))) {

          var nonMutualGenes = 0.0
          for (other <- genes) {
            if (genePatientMatrix.contains(PolimorphismKey(other, sample))) {
              nonMutualGenes += 1
            }
          }

          mutualExScore += 1.0 / nonMutualGenes
        }
      }

      val rank = genes.indexOf(gene).toDouble
      (mutualExScore, mutationCounts, rank)
    })).toList.sortBy(_._2._3)

    val toPrint = new HashMap[String, List[String]]

    val out = new FileWriter(outputPrefix + "_pattern.tsv")

    val header = {
      out.append("Gene")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._1.name))
      out.append("\n").append("Rank")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._2._3.toString()))
      out.append("\n").append("GeneName")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._1.name))
      out.append("\n").append("MutualExclScore")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._2._1.toString()))
      out.append("\n").append("MutQty")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._2._2.toString()))
    }

    val samplesStr = all_samples.map(sample => (sample, {
      val sb = new StringBuilder
      geneMutScoreCounts.foreach(gMsc => sb.append("\t").append(genePatientMatrix.getOrElse(PolimorphismKey(gMsc._1, sample), Polimorphism("Nothing", score = 9)).score.toShort))
      sb.toString
    })).toList.sortBy(_._2)

    val mutualExclusivityPattern = {
      samplesStr.foreach(s => out.append("\n").append(s._1).append(s._2.replaceAll("9", "")))
    }
  }
  
}