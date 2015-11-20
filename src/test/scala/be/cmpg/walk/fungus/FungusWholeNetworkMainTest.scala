package be.cmpg.walk.fungus

import be.cmpg.graph.NetworkReader
import be.cmpg.utils.StatUtils
import scala.collection.Set
import java.nio.file.Paths
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Network
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import java.util.concurrent.Callable
import scala.util.Random
import be.cmpg.walk.Ant
import be.cmpg.graph.interaction.InteractionCostNetworkManager
import java.io.FileWriter
import scala.io.Source
import scala.collection.mutable.HashMap
import java.util.LinkedList
import scala.collection.mutable.ArrayBuffer
import be.cmpg.walk.Path
import collection.immutable.ListMap
import java.io.FileWriter
import be.cmpg.expression.ExpressionNetworkManager

object FungusWholeNetworkMainTest {

  val REAL_SUBNETWORK_SIZE = 5
  val NUMBER_OF_STEPS = 1000
  val FUNGUS_MAX_SIZE = 15

  val selectedGenesSet = new FileWriter(raw"C:\Users\Bram\Desktop/selected.txt")
  val goldenSet = new FileWriter(raw"C:\Users\Bram\Desktop/goldenSet.txt")
  val PossibleGenes = new FileWriter(raw"C:\Users\Bram\Desktop/possibleGenes.txt")

  val interactions = NetworkReader.fromFile(Paths.get("src/test/resources/be/cmpg/graph/network_small_connected.txt"))
  val networkManager = new ExpressionNetworkManager(network = new Network(interactions), evaporation = 0.996)

  new Network(interactions).genes.foreach(gene => PossibleGenes.write(gene + "\n"))
  PossibleGenes.close

  def generateFungus(network: NetworkManager[Gene], score: Double) = {
    (0 to 0).map(_ => networkManager.getGenes.map(gene => {
      new Fungus(
        startGene = gene,
        endGenes = Set(),

        _geneNumberVariable = score,

        network = network)
    })).flatten
  }

  def simulateGeneScores(): (Set[Gene], NodeCostNetworkManager) = {
    /*
     * Select randomly the genes that are going to be important
     */

    var randomImportantGene = networkManager.getGenes.toList(Random.nextInt(networkManager.getGenes.size))
    //val randomImportantGene = Gene("STM4175")
    var resetcounter = 0
    var path = new Path(randomImportantGene)
    while (path.size < REAL_SUBNETWORK_SIZE) {
      if (Random.nextDouble <= 0.2) {
        path.reset
      } else {
        val possibleInteractions = networkManager.getOutgoingInteractionsFor(path.currentEndpoint).filter(interaction => path.canTakeInteraction(interaction)).toList
        if (!possibleInteractions.isEmpty) {
          val next = possibleInteractions(Random.nextInt(possibleInteractions.size))
          path.expand(next)
        } else {
          resetcounter = resetcounter + 1
          if (resetcounter > 20) {
            randomImportantGene = networkManager.getGenes.toList(Random.nextInt(networkManager.getGenes.size))
            path = new Path(randomImportantGene)
            resetcounter = 0
          }
          path.reset
        }
      }
    }

    val importantGenes = path.visitedGenes.toSet

    /*
     * Give scores to all genes.
     */
    val importantGeneMean = 2
    val nonImportantGeneMean = 5

    networkManager.getNetwork().getNodes.foreach(node => node.score = StatUtils.getRandomPoisson(if (importantGenes contains node.gene) importantGeneMean else nonImportantGeneMean))

    importantGenes.foreach(g => println(g + "\t" + networkManager.getNetwork().getNode(g).score))

    println("random scores generated")

    (importantGenes, networkManager)
  }

  def loadGeneScores(): Set[Gene] = {

    for (line <- Source.fromURI(Paths.get("src/test/resources/be/cmpg/graph/simulated_gene_scores.txt").toUri()).getLines()) {
      val splitted = line.split("\t")
      networkManager.getNetwork().getNode(Gene(splitted(0))).score = splitted(1).toDouble
    }

    val importantGenes = Source.fromURI(Paths.get("src/test/resources/be/cmpg/graph/simulated_important_genes.txt").toUri())
      .getLines
      .map(line => Gene(line.split("\t")(0)))
      .toSet

    importantGenes
  }

  def main(args: Array[String]) {
    val importantGenes = simulateGeneScores
    //val importantGenes = loadGeneScores

    /*
     * Start finding the important genes
     */

    val to_json = new HashMap[Gene, ArrayBuffer[Double]]

    val threadpool = java.util.concurrent.Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors())
    //val threadpool = java.util.concurrent.Executors.newSingleThreadExecutor()

    for (i <- 0 to NUMBER_OF_STEPS) {

      val callables = generateFungus(networkManager, FUNGUS_MAX_SIZE).map(walker => {
        new Callable[Option[(Set[Interaction], Double)]] {
          override def call(): Option[(Set[Interaction], Double)] = {
            val subnetwork = walker.selectSubNetwork()
            if (subnetwork.isDefined)
              Some((subnetwork.get, networkManager.scoreSubnetwork(subnetwork.get)))
            else
              return None
          }
        }
      })

      val previous = System.nanoTime()

      val futures = callables.map(threadpool.submit(_))

      val after = System.nanoTime()

      //futures.map(x => println(x.get))
      val paths = futures.map(_.get.get)

      networkManager.updateScores(paths)

      if (i % (NUMBER_OF_STEPS / 50) == 0) {
        print(".")
        for (gene <- importantGenes._1) {
          var values = to_json.get(gene)
          if (values.isEmpty) {
            values = Some(new ArrayBuffer[Double])
            to_json.put(gene, values.get)
          }
          values.get += networkManager.getNetwork().getNode(gene).posteriorProbability
        }

      }

      if (i != NUMBER_OF_STEPS)
        networkManager.evaporate()

    }
    println()
    println("***********************")
    importantGenes._1.foreach(selected => goldenSet.write(selected + "\n"))
    goldenSet.close

    for (gene <- importantGenes._1) {
      var text = "{name:\"" + gene.name + "\",values:["
      var first = true;
      for (value <- to_json(gene)) {
        if (!first)
          text = text + ","
        text = text + "%.3f".format(value * 10)
        first = false
      }
      text = text + "]},"

      //    	println(text)
      /*
    var _TP = 0
    var _FP = 0
    var _TN = 0
    var _FN = 0

    for (gene <- networkManager.getGenes) {
      if (networkManager.probabilityMap(gene) > 0.999) { // positive
        if (importantGenes contains gene) // true positive
          _TP += 1
        else // false positive
          _FP += 1

      } else { // negative
        if (!importantGenes.contains(gene)) // true negative
          _TN += 1
        else
          _FN += 1
      }
    }

    println("***********************")
    println(_TP + "\t" + _FP + "\t" + _TN + "\t" + _FN);
    println("***********************")

    for (gene <- importantGenes) {
      var text = "{name:\"" + gene.name + "\",values:["
      var first = true;
      for (value <- to_json(gene)) {
        if (!first)
          text = text + ","
        text = text + "%.3f".format(value * 10)
        first = false
      }
      text = text + "]},"

      println(text)
      * 
      */
    }

    println("***********************")
    /*
     * Print important gene scores and all gene score evolution during time
     */

    for (gene <- to_json.keys) {
      var text = "{name:\"" + gene.name + "\",values:["
      var first = true;
      for (value <- to_json(gene)) {
        if (!first)
          text = text + ","
        text = text + "%.3f".format(value)
        first = false
      }
      text = text + "]},"

      //      	println(text)
    }

    val resultMap = to_json.map(input => (input._1, input._2.last))
    val orderedMap = ListMap(resultMap.toList.sortBy { _._2 }: _*)
    val selectedGenes = orderedMap.map(input => input._1 + "\t" + input._2.toString).toList.drop(orderedMap.size - 7)
    selectedGenes.foreach(selected => selectedGenesSet.write(selected.split("\t")(0) + "\n"))
    selectedGenesSet.close

    threadpool.shutdown()

    println()
    //networkManager.printResults(System.out)

  }

}