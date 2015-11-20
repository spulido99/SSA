package be.cmpg.expression.time

import be.cmpg.graph.NetworkReader
import be.cmpg.graph.Network
import java.nio.file.Paths
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.walk.fungus.Fungus
import be.cmpg.graph.Gene
import be.cmpg.utils.StatUtils
import java.util.concurrent.Callable
import be.cmpg.graph.Interaction
import scala.collection.Set
import scala.collection.Map
import scala.collection.JavaConversions._
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.mutable.HashMap
import scala.collection.mutable.ListBuffer

object TimeExpressionAnalysis extends App {

  val NUMBER_OF_STEPS = 1000

  val interactions = NetworkReader.fromFile(Paths.get("src/main/resources/be/cmpg/positiveSelectionInBacteria/E.coli K-12 network V1.0.txt"))

  print("Reading expression data. ")
  val timeExpressionMatrix = new CSVReader(new FileReader("src/main/resources/be/cmpg/time/GSE59050.txt"), '\t')
    .readAll()
    .drop(2) // remove first 2 lines
    .map(fields => {
      val gene = fields(1)
      val parentExp = fields(3).toDouble
      val matrix = fields.drop(3).sliding(5, 5).map(_.map(exp => (exp.toDouble - parentExp).abs)).toArray
      (gene, matrix)
    }).toMap

  println("Done.")

  //val threadpool = java.util.concurrent.Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1)
  val threadpool = java.util.concurrent.Executors.newSingleThreadExecutor()

  val network = new Network(interactions)

  val allGenes = network.genes.map(_.name).toSet

  /*
  val networkManager = new TimeExpressionNM(
    network = new Network(interactions),
    timeExpressionMatrix = timeExpressionMatrix,
    all_samples = null,
    mAS_perGene = 1,
    pheromone=0.005,
    evaporation=0.999)
   */

  val toRegress = new CSVReader(new FileReader("src/main/resources/be/cmpg/time/phenotype_regression_tmp.txt"), '\t')
    .readAll()
    .drop(1) // remove first line (header)
    .map(_.drop(1).map(_.toDouble)) // remove the first column (strain name)
    .toArray

  val networkManager = new RegressionTimeExpressionNM(
    network = new Network(interactions),
    timeExpressionMatrix = timeExpressionMatrix,
    toRegress = toRegress,
    all_samples = null,
    mAS_perGene = 1,
    pheromone = 0.005,
    evaporation = 0.996,
    ranked = true)

  for (i <- 0 to NUMBER_OF_STEPS) {

    val walkers = new ListBuffer[Fungus]

    networkManager.getGenes.foreach(gene => {

      val data = timeExpressionMatrix.get(gene.name)

      // create fungus on those genes that in at least one timepoint the change was bigger than 0.2
      if (data.isDefined && data.get.flatten.max > 0.2) {
        walkers += new Fungus(
          startGene = gene,
          endGenes = Set(),
          _geneNumberVariable = 3,//2 + StatUtils.getRandomPoisson(2), // ?? Random size or static 3 size????
          network = networkManager)
      }
    })

    val callables = walkers.map(walker => {
      new Callable[Option[(Set[Interaction], Double)]] {
        override def call(): Option[(Set[Interaction], Double)] = {
          val subnetwork = walker.selectSubNetwork()
          val score = networkManager.scoreSubnetwork(subnetwork.get)
          Some((subnetwork.get, score))
        }
      }
    })

    val previous = System.nanoTime()

    val futures = callables.map(threadpool.submit(_))

    val after = System.nanoTime()

    //futures.map(x => println(x.get))
    val paths = futures.map(_.get.get)

    networkManager.updateScores(paths)

    if (i % (NUMBER_OF_STEPS / 50).ceil == 0) {
      print(".")
    }

    if (i != NUMBER_OF_STEPS)
      networkManager.evaporate()

  }

  threadpool.shutdown()

  val rankedGenes = networkManager.getRankedAllGenes
  var count = 0
  println("*******************")
  for (gene <- rankedGenes) {
    //if (networkManager.getPosteriorProbability(gene) > 0.99) {
    count += 1
    println(gene.name + "\t" + count + "\t" + networkManager.getPosteriorProbability(gene))
  }
  println()
  println("-----\t" + count)
}