package be.cmpg.cluster

import java.io.FileWriter
import org.json.JSONObject
import org.json.JSONArray
import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.NodeCostNetworkManager
import scala.collection.JavaConversions._
import scala.io.Source
import java.io.FileInputStream

object PValuesPrintPattern extends App {

  val parser = ArgumentsParser.getBasicArgParser("SSA.ME Print patterns from a list of genes")

  parser.opt[String]("geneList") action { (x, c) =>
    c.copy(genes = x.split(",").toList)
  } text ("Gene list to print network and pattern (comma separated).")

  val configOpt = parser.parse(args, PvalueConfig(

    
     //Default Values
     
    iterations = 1000,
    reinforcement = 0.005,
    forgetfulness = 0.996))

  if (configOpt.isEmpty) {
    parser.showUsageAsError
    System.exit(1)
  }
  // parser.parse returns Option[C]

  val config = configOpt.get

  val (interactions, network, pValueMap, geneList) = DataLoader.loadData(config)

  val networkManager = new ClusterNetworkManager(
    network = network,
    genePValueMatrix = pValueMap,
    evaporation = 0,
    pheromone = 0,
    ranked = false)

  printPattern(config.outputPrefix + ".selected", config.genes.map { n => Gene(n) }, networkManager, pValueMap, config.outputFolder,config.inputFolder)

  def printPattern(outputPrefix: String, genes: List[Gene], networkManager: NodeCostNetworkManager, pValueMap: Map[String, Double], outputFolder: String,inputFolder:String) = {

    val results = new JSONObject

    val edges = new JSONArray

    val out_edges = new FileWriter(outputFolder + outputPrefix + "_edges.tsv")
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
            .accumulate("evidence", i.evidence.mkString(","))
          edges.put(edgeInfo)
          out_edges.write(i.from.name + "\t" + i.to.name + "\t" + score + "\n")
        }
      }

    out_edges.close()

    val nodes = new JSONArray

    val out_nodes = new FileWriter(outputFolder + outputPrefix + "_nodes.tsv")
    out_nodes.write("GeneSymbol\tPosteriorP\tSelected\tBestSSN\n")
    genes
      .foreach(g => {

        println(g.name)

        val geneInfo = new JSONObject()
          .accumulate("name", g.name)
          .accumulate("selected", genes.contains(g))

        nodes.put(geneInfo)
        try {
          val node = networkManager.getNetwork().getNode(g)
          out_nodes.write(g.name + "\t" + node.posteriorProbability + "\t" + genes.contains(g) + "\t" + node.convergenceIteration + "\t"  + node.bestSubnetwork._1 + "\t" + node.bestSubnetwork._2 + "\n")
        } catch {
          case t: Throwable => out_nodes.write(g.name + "\t" + "NA" + "\t" + genes.contains(g) + "\t" + "NA" +  "\t" + "NA" + "\n")
        }
      })

    out_nodes.close()

    results.put("edges", edges)
    results.put("nodes", nodes)

    val html = Source.fromInputStream(new FileInputStream(inputFolder + "adjustedNetworkOutput.html")).mkString

    val out_json = new FileWriter(outputFolder + outputPrefix + "_network.html")
    out_json.write(html.replace("$${graph}$$", results.toString.replace("\"[", "[").replace("]\"", "]")))

    out_json.close()

  }

}