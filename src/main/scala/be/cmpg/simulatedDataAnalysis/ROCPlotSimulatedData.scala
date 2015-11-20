package be.cmpg.simulatedDataAnalysis

import java.io.FileWriter
import be.cmpg.walk.Path
import be.cmpg.graph.interaction.NodeCostNetworkManager
import be.cmpg.graph.Gene
import scala.util.Random
import org.apache.commons.math3.distribution.ChiSquaredDistribution
import be.cmpg.graph.NetworkReader
import java.nio.file.Paths
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import java.util.concurrent.Callable
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.walk.fungus.Fungus
import scala.collection.Set
import java.util.concurrent.Executor
import java.util.concurrent.ExecutorService
import scala.collection.mutable.HashSet
import org.apache.commons.math3.distribution.BetaDistribution
import scala.collection.mutable.ListBuffer
import scala.collection.mutable.HashMap
import be.cmpg.expression.ExpressionNetworkManager
import be.cmpg.expression.ExpressionNetworkManager
import org.apache.commons.math3.distribution.NormalDistribution

object ROCPlotSimulatedData extends App {

  val repeats = 100
  
  def simulateGeneScores(networkManager: NodeCostNetworkManager, subnetworkSize: Int, chiDfRest: Double = 5): (Set[Gene], HashMap[Gene, Double]) = {
    /*
     * Select randomly the genes that are going to be important
     */

    var randomImportantGene = networkManager.getGenes.toList(Random.nextInt(networkManager.getGenes.size))
    //val randomImportantGene = Gene("STM4175")
    var resetcounter = 0
    var path = new Path(randomImportantGene)
    while (path.size < subnetworkSize) {
      if (Random.nextDouble <= 0.2) {
        path.restart
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
    val importantDist: BetaDistribution = new BetaDistribution(2, 1.5) // Important
    val restDist: NormalDistribution = new NormalDistribution(0, 0.3)

    val geneExpression = new HashMap[Gene, Double]
    
    networkManager.getGenes.foreach(g => geneExpression.put(g,
      if (importantGenes contains g)
        importantDist.inverseCumulativeProbability(Random.nextDouble)
      else
        math.min(math.abs(restDist.inverseCumulativeProbability(Random.nextDouble)), 1)))

    importantGenes.foreach(g => println(g + "\t" + geneExpression(g)))

    println("random scores generated")

    (importantGenes, geneExpression)
  }

  def generateFungus(networkManager: NodeCostNetworkManager, score: Double) = {
    
    val walkers = networkManager.getNetwork().getNodes()
      .filter { _.score > 0.5 }
      .map { node => {
        new Fungus(
                startGene = node.gene,
                endGenes = Set(),
                _geneNumberVariable = 5,
                network = networkManager)
      } }
    
    //walkers.view.toList
    walkers
  }

  def getRankedGenesByPhac(networkManager: ExpressionNetworkManager, numberOfSteps: Int = 1000, maxFungusSize: Double = 15, threadpool: ExecutorService) = {
    for (i <- 0 to numberOfSteps) {

      val callables = generateFungus(networkManager, maxFungusSize).map(walker => {
        new Callable[Option[(Set[Interaction], Double)]] {
          override def call(): Option[(Set[Interaction], Double)] = {
            val subnetwork = walker.selectSubNetwork()
            if (subnetwork.isDefined) {
            	Some((subnetwork.get, networkManager.scoreSubnetwork(subnetwork.get)))
            } else
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

      if (i % 50 == 0) print(".")

      if (i != numberOfSteps)
        networkManager.evaporate()

    }
    println()

    networkManager.getRankedAllGenes

  }

  val threadpool = java.util.concurrent.Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors() - 1)
  //val threadpool = java.util.concurrent.Executors.newSingleThreadExecutor()

  val networkFiles = Map(
    "small" -> List("src/test/resources/be/cmpg/graph/network_small_connected.txt"),
    "bacteria" -> List("src/main/resources/be/cmpg/graph/network_1.txt"),
    "human" -> List("src/main/resources/be/cmpg/graph/human/edges.csv", "src/main/resources/be/cmpg/graph/human/nodes.csv")
	);

  for (networkName <- networkFiles.keys) {
    val interactions = if (networkFiles(networkName).size == 1)
				       	   NetworkReader.fromFile(Paths.get(networkFiles(networkName)(0)))
				       else
				    	   NetworkReader.fromCytoScapeFiles(Paths.get(networkFiles(networkName)(0)), Paths.get(networkFiles(networkName)(1)))._1

    val network = new Network(interactions)

    val allGenes = network.genes.map(_.name).toSet
    val size = 10

    val writer = new FileWriter(networkName+size+"_output.txt")
    writer.write("Rep\tN\tphac_aspecificity\tphac_sensitivity\tneighbourhood_aspecificity\tneighbourhood_sensitivity\tstat_aspecificity\tstat_sensitivity\n")

    for (r <- 0 to repeats) {

      println("---> " + r)

      val networkManager = new ExpressionNetworkManager(network = network)
      
      println("- Network Manager created")
      
      val (importantGenesR, geneExpression) = simulateGeneScores(networkManager, size)
      
      val importantGenes = importantGenesR.map(_.name).toSet
      networkManager.geneExpression = geneExpression.view.toMap;
      
      println("- Important genes simulated")
      /*
	   * Ranked genes by score
	   */
      print("- Rank by Score. ")
      val rankedGenesByScore = networkManager.getGenes.toList.sortBy(g => networkManager.geneExpression(g))(Ordering[Double].reverse).toList.map(_.name)
      println("Done!")

      /*
	   * Ranked genes by Neighbourhood
	   */
      print("- Rank by neighbourhood. ")
      val neighbourhoodRank = new NeighbourhoodScoringBySeedNode(networkManager).rankedGenes.reverse
      println("Processing. Done!")
      /*
	   * Ranked genes by phac
	   */
      print("- Start rank by Phac. ")
      val rankedGeneByPhac = getRankedGenesByPhac(networkManager, threadpool = threadpool).toList.map(_.name);
      println("Done!")

      val phacSelected = new HashSet[String]
      val neighbourhoodSelected = new HashSet[String]
      val statSelected = new HashSet[String]
      
      print("calculate ROC values. ")
      for (n <- 0 to allGenes.size - 1) {

        phacSelected += rankedGeneByPhac(n)
        val phacROCValues = new ROCPlotValuesCreator(phacSelected.view.toSet, importantGenes, allGenes).get

        neighbourhoodSelected += neighbourhoodRank(n)
        val neighbourhoodROCValues = new ROCPlotValuesCreator(neighbourhoodSelected.view.toSet, importantGenes, allGenes).get

        statSelected += rankedGenesByScore(n)
        val statROCValues = new ROCPlotValuesCreator(statSelected.view.toSet, importantGenes, allGenes).get

        writer.write(r + "\t" + n + "\t" + phacROCValues.sensitivity + "\t" + phacROCValues.aspecificity + "\t" + neighbourhoodROCValues.sensitivity + "\t" + neighbourhoodROCValues.aspecificity + "\t" + statROCValues.sensitivity + "\t" + statROCValues.aspecificity + "\n")
      }
      val rankedGenes = networkManager.getRankedAllGenes
      var count = 0
      for (gene <- rankedGenes) {
        val score = networkManager.getPosteriorProbability(gene);
        if (score > 0.99)
          count += 1
      }
      println()
      println("-----\tRanked Genes:" + count)
      
      println("Done!")
      
      
      
    }

    writer.close()

  }

  threadpool.shutdown()

}