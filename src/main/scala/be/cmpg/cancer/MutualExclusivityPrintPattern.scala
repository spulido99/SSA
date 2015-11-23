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

object MutualExclusivityPrintPattern extends App {

  case class Config(
    outputPrefix: String,
    refNetwork: Seq[String],
    genes: List[String] = List(),
    maf: Option[File] = None,
    expression: Map[String, File] = Map())

  val parser = new scopt.OptionParser[Config]("SSA.ME.") {
    head("Mutual exclusivity analysis by Small Subnetwork Analysis", "0.1")

    opt[Seq[String]]('n', "refNetwork") action { (x, c) =>
      c.copy(refNetwork = x)
    } text ("The reference network (HT, hiII14, ppi, cell-map, nci-nature, reactome or/and files. HT_hiII14 used by default.)")

    opt[String]('o', "outputPrefix") action { (x, c) =>
      c.copy(outputPrefix = x)
    } text ("The prefix for the analisis output files")

    opt[String]('l', "genelist") required () action { (x, c) =>
      c.copy(genes = x.split(",").toList)
    } text ("Gene list to print network and pattern (comma separated).")

    opt[File]('m', "maf") required () action { (x, c) =>
      c.copy(maf = Some(x))
    } text ("Mutation .maf file")

    opt[Map[String, File]]('e', "expression") action { (x, c) =>
      c.copy(expression = x)
    } text ("expression files from GISTIC (... -e cnv_peaks=<file1>,exp=<file2>,cnv_thresholds=<file3> ...).")

    help("help")

  }
  // parser.parse returns Option[C]

  val helper = new CancerHelper

  val configOption = parser.parse(args, Config(
    /*
       * Default Values
       */
    refNetwork = Seq("HT", "hiII14", "reactome"),
    outputPrefix = "ME"))

  if (configOption.isEmpty) {
    parser.showUsageAsError
    System.exit(1)
  }

  val config = configOption.get

  val genePatientMatrix = {

    val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]
    helper.loadMaf(genePatientMatrix, config.maf.get)

    if (config.expression.size == 3)
      helper.loadExpression(genePatientMatrix,
        config.expression("cnv_peaks"),
        config.expression("exp"),
        config.expression("cnv_thresholds"))

    genePatientMatrix.toMap
  }

  val interactions = helper.loadNetwork(config.refNetwork)
  val network = new Network(interactions)

  val networkManager = new MutualExclusivityNetworkManager(
    network = network,
    genePatientMatrix = genePatientMatrix,
    minimumSamplesAltered = 0,
    pheromone = 0,
    evaporation = 0,
    ranked = false)

  printPattern(config.outputPrefix, config.genes.map { n => Gene(n) }, networkManager, genePatientMatrix)

  def printPattern(outputPrefix: String, genes: List[Gene], networkManager: NodeCostNetworkManager, genePatientMatrix: Map[PolimorphismKey, Polimorphism], additional: Map[Gene, Map[String, Any]] = Map()) = {

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

    
    val cgc = new CSVReader(new FileReader("cancergenelists/Census_allFri20150102.csv")) // CGC: Cancer Gene Census
                                .readAll()
                                .drop(1)
                                .map(_(0))
                                .toSet
                                
    val malacard = new CSVReader(new FileReader("cancergenelists/malacardBreastCancer.tsv"), '\t') // Malacards
                                .readAll()
                                .drop(1)
                                .map(_(1))
                                .toSet            
    
   val ncg = new CSVReader(new FileReader("cancergenelists/cancergenes.txt"), '\t') // NCG 4.0 : Network of Cancer Genes
                                .readAll()
                                .drop(1)
                                .map(_(1))
                                .toSet
                                
                                
    def getKnownCancerGene(gene:Gene) = {
      if (malacard.contains(gene.name)) "malacard"
      else if (cgc.contains(gene.name)) "cgc"
      else if (ncg.contains(gene.name)) "ncg"
      else "unkown"
    }  
    
    val truePositiveGeneList = cgc ++ malacard ++ ncg;
    
    val nodes = new JSONArray

    val out_nodes = new FileWriter(outputPrefix + "_nodes.tsv")
    out_nodes.write("GeneSymbol\tPosteriorP\tSelected\tConvergenceIter\tKnownCancerGene\t\tBestSSN\n")
    println("GeneSymbol\tMalacard\tCGC\tNCG")
    genes
      .foreach(g => {
        
        println(g.name+"\t"+malacard.contains(g.name)+"\t"+cgc.contains(g.name)+"\t"+ncg.contains(g.name))
        
        val geneInfo = new JSONObject()
          .accumulate("name", g.name)
          .accumulate("selected", genes.contains(g))
          .accumulate("knownCancerGene", getKnownCancerGene(g))
          
        val additionalInfo = additionalNodeInfo.get(g)
        if (additionalInfo.isDefined) {
          additionalInfo.get.foreach { e =>
            geneInfo.accumulate(e._1, e._2)
          }
        }
        nodes.put(geneInfo)
        try {
        	val node = networkManager.getNetwork().getNode(g)
        	out_nodes.write(g.name + "\t" + node.posteriorProbability+ "\t" + genes.contains(g) + "\t" + node.convergenceIteration + "\t" + truePositiveGeneList.contains(g.name) + "\t" + node.bestSubnetwork._1+"\n")
        } catch {
          case t: Throwable => out_nodes.write(g.name + "\t" + "NA" + "\t" + genes.contains(g) + "\t" + "NA" + "\t" + truePositiveGeneList.contains(g.name) + "\t" + "NA" + "\n")
        }
      })

    out_nodes.close()

    results.put("edges", edges)
    results.put("nodes", nodes)

    val baseNetworkOutput = scala.io.Source.fromFile("baseNetworkOutput.html")
    val lines = baseNetworkOutput.mkString
    baseNetworkOutput.close()

    val out_json = new FileWriter(outputPrefix + "_network.html")
    out_json.write(lines.replace("$${graph}$$", results.toString))

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