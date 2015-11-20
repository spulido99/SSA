package be.cmpg.cancer

import java.io.FileReader
import java.nio.file.Paths
import java.util.concurrent.Callable
import scala.collection.JavaConversions.asScalaBuffer
import scala.collection.JavaConversions.mutableSeqAsJavaList
import scala.collection.Set
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import scala.collection.mutable.LinkedList
import au.com.bytecode.opencsv.CSVReader
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.graph.NetworkReader
import be.cmpg.walk.fungus.Fungus
import be.cmpg.utils.StatUtils
import be.cmpg.graph.Interaction
import java.io.File
import java.io.PrintWriter
import scala.collection.mutable.StringBuilder
import be.cmpg.walk.SubNetworkSelector
import java.io.FileWriter

object MemoComparissonAnalysis extends App {

  if (this.args.size < 1)
    throw new IllegalArgumentException("Usage: java -jar phac-assembly-0.1.jar <number of iterations:100> <pheromones:0.01> <evaporation:0.996> <reference networks:HT_hiII or MEMo> <useRank:true or false> <maxOutputGenes:200 genes for output network>")

  val iterations = this.args(0).toInt
  val pheromones = this.args(1).toDouble
  val evaporation = this.args(2).toDouble
  val refNetwork = this.args(3)
  val useRank = this.args(4).toBoolean
  val maxOutputGenes = this.args(5).toInt
  val mAS_perGene = 3

  /*val (interactions, translateGenesToEntrez) = NetworkReader.fromCytoScapeFiles(
    Paths.get("src/main/resources/be/cmpg/graph/human/edges.csv"),
    Paths.get("src/main/resources/be/cmpg/graph/human/nodes.csv"))*/

  val helper = new CancerHelper()

  val interactions = helper.loadReferenceNetwork(refNetwork)

  val network = new Network(interactions)

  val translateGenesToEntrez: Map[String, String] = network.genes.map(g => (g.name, g.name)).toMap

  //val genePatientMatrix = helper.loadFromDendrixFiles()
  //val genePatientMatrix = helper.loadFromTcgaFiles()
  val genePatientMatrix = helper.loadFromTcga2012Files()
  val all_samples = genePatientMatrix.map(_._1.sample).toSet
  val geneList = new HashSet[Gene]
  network.genes.foreach(gene => {
    var count = 0
    for (sample <- all_samples) {
      if (genePatientMatrix contains (PolimorphismKey(gene, sample))) {
        count += 1
      }
    }

    if (count > 10) {
      geneList.add(gene)
    }
  })

  /*new CSVReader(new FileReader("src/test/resources/dendrix/GBM.glst"), '\t')
    .readAll()
    .foreach(fields => {
      val hugoGeneName = fields(0)
      val entrezId = Gene(translateGenesToEntrez.get(hugoGeneName).getOrElse(hugoGeneName))
      if (network.genes.contains(entrezId)) {
        geneList.add(entrezId)
      }
    })*/

  val networkManager = new MutualExclusivityNetworkManager(
    network = network,
    genePatientMatrix = this.genePatientMatrix,
    minimumSamplesAltered = mAS_perGene,
    pheromone = pheromones,
    evaporation = evaporation,
    ranked = useRank)

  val walkers: Set[SubNetworkSelector] = geneList.map(gene =>
    //new Fungus(
    new Hi2iiFungus(
      startGene = gene,
      geneNumberVariable = 2 + StatUtils.getRandomPoisson(1),
      //geneNumberVariable = 3 + StatUtils.getRandomPoisson(0.5), // 0 - 60%, 1 - 30%, 2 - 10%
      //geneNumberVariable = 3,
      network = networkManager)).toSet

  networkManager.run(iterations, geneList.view.toSet, Runtime.getRuntime().availableProcessors() / 2, defwalkers = Some(walkers))

  val cosmicGenes = new CSVReader(new FileReader(Paths.get("src/test/resources/cosmic/Census_allFri20150102.csv").toFile()))
    .readAll()
    .drop(1)
    .map(_(2))
    .toSet;

  val rankedGenes = networkManager.getRankedAllGenes.take(maxOutputGenes).toSet
  // MEMo Ranked Genes

  /*val rankedGenes = Set(
      Gene("AKT1"), Gene("AKT3"), Gene("ATM"), Gene("BRCA1"), 
      Gene("CCND1"), Gene("CDKN2A"), Gene("CHEK2"), Gene("EGFR"), 
      Gene("ERBB2"), Gene("IGF1R"), Gene("IKBKB"), Gene("MAP2K4"), 
      Gene("MAP3K1"), Gene("MDM2"), Gene("MDM4"), Gene("MYC"), 
      Gene("NBN"), Gene("PAK1"), Gene("PIK3CA"), Gene("PIK3R1"), 
      Gene("PTEN"), Gene("RB1"), Gene("TP53"))
*/
  val translateEntrezToGene = translateGenesToEntrez.map(x => (x._2, x._1)).toMap
  println()
  println("********* Selected Subnetwork size: " + rankedGenes.size)

  val geneMutScoreCounts = rankedGenes.map(gene => (gene, {
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
        for (other <- rankedGenes) {
          if (genePatientMatrix.contains(PolimorphismKey(other, sample))) {
            nonMutualGenes += 1
          }
        }

        mutualExScore += 1.0 / nonMutualGenes
      }
    }

    val rank = networkManager.rankedGenes.view.toList.indexOf(gene).toDouble
    (mutualExScore, mutationCounts, rank)
  })).toList.sortBy(_._2._3)

  val out_edges = new FileWriter("edges."+refNetwork+".tsv")

  networkManager.getAllInteractions
    .filter { i => i.to != i.from && rankedGenes.contains(i.from) && rankedGenes.contains(i.to) }
    .foreach { i => out_edges.write(i.from.name + "\t" + i.to.name + "\t" + networkManager.scoreSubnetwork(Set(i)) + "\n") }

  out_edges.close()

  val out_nodes = new FileWriter("nodes."+refNetwork+".tsv")
  
  rankedGenes
    .zipWithIndex
    .foreach( gi => out_nodes.write(gi._1.name + "\t" + gi._2 + "\t" + networkManager.getPosteriorProbability(gi._1) + "\n"))
  
  out_nodes.close()
  
  val toPrint = new HashMap[String, List[String]]

  val header = {
    val sb = new StringBuilder
    sb.append("Gene")
    geneMutScoreCounts.foreach(g => sb.append("\t").append(g._1.name))
    sb.append("\n").append("Rank")
    geneMutScoreCounts.foreach(g => sb.append("\t").append(g._2._3))
    sb.append("\n").append("GeneName")
    geneMutScoreCounts.foreach(g => sb.append("\t").append(translateEntrezToGene.getOrElse(g._1.name, g._1.name)))
    sb.append("\n").append("MutualExclScore")
    geneMutScoreCounts.foreach(g => sb.append("\t").append(g._2._1))
    sb.append("\n").append("MutQty")
    geneMutScoreCounts.foreach(g => sb.append("\t").append(g._2._2))
    sb.append("\n").append("IsCosmic")
    geneMutScoreCounts.foreach(g => sb.append("\t").append(cosmicGenes.contains(g._1.name)))
    sb.toString
  }
  print(header)

  val samplesStr = all_samples.map(sample => (sample, {
    val sb = new StringBuilder
    geneMutScoreCounts.foreach(gMsc => sb.append("\t").append(genePatientMatrix.getOrElse(PolimorphismKey(gMsc._1, sample), Polimorphism("Nothing", score = 9))._type))
    sb.toString
  })).toList.sortBy(_._2)

  val mutualExclusivityPattern = {
    val sb = new StringBuilder
    samplesStr.foreach(s => sb.append("\n").append(s._1).append(s._2.replaceAll("9", "")))
    sb.toString
  }
  print(mutualExclusivityPattern)

  /*
  println("*******************")
  networkManager.acumulativeInteractionsScore.toList
    .map(v => {
      /*var text = v._1.toString
  	  text += "\t" + v._1.genes.map(networkManager.getGeneScore(_)).min
  	  text += "\t" + v._2._1
  	  text += "\t" + v._2._2
  	  text += "\t" + v._2._2
  	  text += "\t" + v._2._1 / v._2._2
  	    */
      (v._1, v._2._1 / v._2._2)
      //(v._1, v._1.genes.map(networkManager.getGeneScore(_)).min * v._2._1 / v._2._2)
      //(v._1, v._1.genes.map(networkManager.getGeneScore(_)).min * v._2._2)
    })
    .foreach(e => if (e._2 > 0.0) println(e._1 + "\t" + e._2))
    * 
    */

}