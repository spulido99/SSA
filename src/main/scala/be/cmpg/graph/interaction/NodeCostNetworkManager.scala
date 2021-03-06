package be.cmpg.graph.interaction

import scala.collection.mutable.HashMap
import scala.util.Random
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.walk.SubNetworkSelector
import scala.collection.Traversable
import scala.collection.mutable.ArrayBuffer
import scala.collection.Map
import scala.collection.Set
import scala.collection.mutable.LinkedHashSet
import be.cmpg.graph.Interaction
import scala.collection.mutable.LinkedHashMap
import scala.Option
import be.cmpg.graph.Interaction
import be.cmpg.walk.fungus.Fungus
import java.util.concurrent.Callable
import be.cmpg.utils.StatUtils
import java.util.LinkedList
import scala.collection.mutable.HashSet
import scala.util.control.Breaks

abstract class NodeCostNetworkManager(network: Network,
                                      pheromone: Double = 0.005,
                                      evaporation: Double = 0.996,
                                      ranked: Boolean = false,
                                      minProb: Double = 0.01,
                                      initialProb: Double = 0.50,
                                      convergenceThreshold: Double = 0.99) extends NetworkManager[Gene](network) {

  def getMinProb = minProb

  //val acumulativeInteractionsScore = new HashMap[Interaction, (Double, Integer)] // interaction -> (sum of scaled scores, number of times selected)
  val rankedGenes = new LinkedHashSet[Gene]
  //val bestSubnetworkRankedGenes = new HashMap[Gene, (Set[Interaction], Double)]
  //val probabilityMap = initProbabilityMap()
  //var probabilityHistory:Map[Gene, LinkedList[Double]] = _

  /*
   * Initialize node probabilities
   */
  network.getNodes().foreach { node => node.posteriorProbability = initialProb }

  def getRankedGenes() {

  }

  def getRankedAllGenes() = {
    val toReturn = new LinkedHashSet[Gene]

    toReturn ++= rankedGenes

    network.getNodes()
      .toList
      .sortBy(node => -node.posteriorProbability)
      .foreach { node => if (!rankedGenes.contains(node.gene)) toReturn += node.gene }

    toReturn.view
  }

  /*
   * Score subnetwork by the number of times a gene exist in subnetworks
   * 
  def scoreSubnetwork(subnetworks: Traversable[Set[Interaction]]) = {
    val genes = HashMap[Gene, Double]()
    for (subnetwork <- subnetworks) {
      for (interaction <- subnetwork) {
        if (!genes.contains(interaction.from)) genes.put(interaction.from, 1.0) else genes.put(interaction.from, genes(interaction.from) + 1.0)
        if (!genes.contains(interaction.to)) genes.put(interaction.to, 1.0) else genes.put(interaction.to, genes(interaction.to) + 1.0)
      }
    }
    genes
  }
  *
  */

  def debugForGenes(genes: Set[Gene], iterations: Int = 1) = {
    debug = Some(genes, iterations)
  }

  var debug: Option[(Set[Gene], Int)] = None
  var iterationsLeft = 1000;

  private var fakeMax: Option[Double] = None
  /**
   * Use the stated max instead of the actual max.
   *
   * Useful when the max of the subnetwork scoring function is know
   */
  def setUserFakeMax(fakeMax: Double) = this.fakeMax = Some(fakeMax)

  def scaleSubnetworksScoresByGene(inputSubnetworkScores: Traversable[WalkerResult]): Map[Gene, Double] = {

    
    if (debug.isDefined) {
      
      val results = inputSubnetworkScores.filter { result => debug.get._1 contains result.walker.getStartGene() }
                                          .groupBy { _.walker.getStartGene() }
                                          .mapValues { _.maxBy { _.score } }
      
      if (iterationsLeft == 1000) {
        println()
        debug.get._1.foreach(g => print("\tNetScore" + g.name))
        debug.get._1.foreach(g => print("\tNodeProb" + g.name))
        debug.get._1.foreach(g => print("\tSelectedNet" + g.name))
      }
  
      if (debug.isDefined && iterationsLeft % debug.get._2 == 0) {
  
        println()
        print(iterationsLeft)
        debug.get._1.foreach(g => print("\t" + results.getOrElse(g, WalkerResult(null, Set(), 0.0)).score))
        debug.get._1.foreach(g => print("\t" + network.getNode(g).posteriorProbability))
        debug.get._1.foreach(g => print("\t" + results.getOrElse(g, WalkerResult(null, Set(), 0.0)).subnetwork.map(_.genes).flatten.map(_.name).toSet))
  
      }
    }
    
    /*
     * filter out NaN values
     */
    var subnetworkScores = inputSubnetworkScores.filter(sn => !sn.subnetwork.isEmpty && !sn.score.isNaN())

    if (ranked) {
      var rank = 0.0;
      subnetworkScores = subnetworkScores
        .toList
        .sortBy(_.score)
        .filter(_.score > 0)
        .map(s => {
          rank += 1
          WalkerResult(s.walker, s.subnetwork, rank)
        })
    }

    var minScore = Double.MaxValue
    var maxScore = Double.MinValue

    /*
     * store score by gene
     */
    val scores = HashMap[Gene, WalkerResult]()
    //val scores_i = HashMap[Interaction, Double]()
    for (result <- subnetworkScores) {

      /*
       * update the min/max scores
       */
      minScore = math.min(minScore, result.score)
      maxScore = math.max(maxScore, result.score)

      /*
       * Set to each gene the highest score found in all subnetworks it is present
       */
      for (interaction <- result.subnetwork) {
        /*
         *  add If the gene (from) is not present, 
         *  or the stored score is lower than the current score
         */
        val fromNode = network.getNode(interaction.from)

        val cScoreFrom = scores.get(interaction.from)
        if (!cScoreFrom.isDefined || result.score > cScoreFrom.get.score) {
          scores.put(interaction.from, result)

          if (fromNode.bestSubnetwork._2 < result.score) {
            fromNode.bestSubnetwork = (result.subnetwork, result.score)
          }
        }

        /*
         *  add If the gene (to) is not present, 
         *  or the stored score is lower than the current score
         */
        val cScoreTo = scores.get(interaction.to)
        val toNode = network.getNode(interaction.to)
        if (!cScoreTo.isDefined || result.score > cScoreTo.get.score) {
          scores.put(interaction.to, result)

          if (toNode.bestSubnetwork._2 < result.score) {
            toNode.bestSubnetwork = (result.subnetwork, result.score)
          }
        }

        /*
         * 
        
        val cScoreI = scores_i.get(interaction)
        if (!cScoreI.isDefined || score > cScoreI.get) {
          scores_i.put(interaction, score)
        }
         */
      }
    }

    if (fakeMax.isDefined) {
      if (fakeMax.get < maxScore)
        throw new RuntimeException("Actual maximum score is greater than fake max score")
      maxScore = fakeMax.get
    }

    /*
     * Scale scores between min/max
     */
    val range = maxScore - minScore
    val scaledScores = scores.map(x => {
      if (range == 0) {
        (x._1, 0.0)
      } else {
        val scaledValue = (x._2.score - minScore) / (range);
        (x._1, scaledValue)
      }
    }).toMap

    iterationsLeft -= 1

    scaledScores
  }

  /*
   * (subnetwork, score)
   */
  override def updateScores(subnetworkScores: Traversable[WalkerResult]): Map[Gene, Double] = {

    /*
     * Update by the best scored subnetwork in which the gene is part of. 
     */
    val genesScores = scaleSubnetworksScoresByGene(subnetworkScores)

    genesScores.foreach(g => {

      val node = network.getNode(g._1)

      val increase = 1 + pheromone * g._2
      val newValue = node.posteriorProbability * increase

      if (!rankedGenes.contains(g._1) && newValue >= 0.95) {
        rankedGenes.add(g._1)
      }

      node.posteriorProbability = math.min(1, math.max(newValue, minProb)) // Do not allow any interation to go over 1 or less than 0.01
    })

    genesScores
  }

  def evaporate(geneScores: Map[Gene, Double] = Map()) = {

    /*
     * Given that all genes belonged to at least one path, re-center them at 0.5
     */

    //var sum = 0.0
    var min = Double.MaxValue
    var max = Double.MinValue
    network.getNodes().foreach(node => {
      node.posteriorProbability = ((node.posteriorProbability * evaporation) min 1.0) // check to not minProv <= value <= 1.0
      //sum += node.posteriorProbability
      min = min.min(node.posteriorProbability)
      max = max.max(node.posteriorProbability)
    })
    /*
     * Given that all genes belonged to at least one path, re-center them at 0.5
     */
    val mid = (min + max) / 2.0
    //val mid = sum / network.getNodes().size

    network.getNodes().foreach(node => {
      node.posteriorProbability = (((node.posteriorProbability - mid + initialProb) min 1.0))
      //node.posteriorProbability = ((node.posteriorProbability min 1.0) max minProb)
    })
  }

  override def getPosteriorProbability(gene: Gene): Double = {
    network.getNode(gene).posteriorProbability
  }

  def converged(): Boolean = { // iteration:Int

    var convergedNodes = network.getNodes().filter { _.convergenceIteration > 0 }

    /*if (iteration % 50.ceil == 0) {
       println(iteration + " >> "+ convergedNodes.size.toDouble + "/" + network.getNodes.size.toDouble + " = "+convergedNodes.size.toDouble / network.getNodes.size.toDouble)
    }*/

    convergedNodes.size.toDouble / network.getNodes.size.toDouble > convergenceThreshold
  }

  override def getRandomInteraction(selector: SubNetworkSelector): Option[Interaction] = {
    val listChildInteractions = selector.getPossibleInteractions()
    if (listChildInteractions.isEmpty) {
      /*
       * there are no more posible child interactions.
       * Return max value to finish loop calling function
       */
      None
    } else {

      val posteriorProbs = new ArrayBuffer[Double](listChildInteractions.size)
      for (i <- listChildInteractions) {
        val targetGene = if (selector.getVisitedGenes.contains(i.from)) i.to else i.from
        val posteriorProb = math.max(getPosteriorProbability(targetGene), minProb)
        posteriorProbs += posteriorProb
      }

      var ran = Random.nextDouble * posteriorProbs.sum
      var i = 0

      val size = listChildInteractions.size
      while (ran > 0 && size > i) {
        ran -= math.max(posteriorProbs(i), minProb)
        i += 1
      }

      Some(listChildInteractions(i - 1))
    }
  }

  def run(iterations: Int, startingNodes: Set[Gene], processors: Int = 1, endNodes: Set[Gene] = Set(), defwalkers: Option[Set[SubNetworkSelector]] = None, storeNodePHistory: Boolean = false, storeNodePHistIter: Int = 50, subnetworkSize: Option[Int] = None) = {

    iterationsLeft = iterations

    if (storeNodePHistory) {
      network.getNodes.foreach(_.probabilityHistory = new LinkedList[Double])
    }

    val threadpool = if (debug.isDefined || processors == 1)
      java.util.concurrent.Executors.newSingleThreadExecutor()
    else
      java.util.concurrent.Executors.newFixedThreadPool(processors)

    val walkers = if (defwalkers.isDefined)
      defwalkers.get
    else
      startingNodes.map(gene =>
        new Fungus(
          startGene = gene,
          endGenes = endNodes,
          network = this)).toSet

    // Initialize convergence iterations. This is not needed for the mutualExclusivityAnalysis but is needed for the bootstrapCalculator as it otherwise stops after any solution has converged.
    network.getNodes().foreach(_.convergenceIteration = -1)

    for (iteration <- 0 to iterations if (!converged())) {

      if (debug.isEmpty && iteration % (iterations / 5).ceil == 0) {
        println()
        print((iterations - iteration) + "\t- ")
      }
      if (debug.isEmpty && iteration % (iterations / storeNodePHistIter).ceil == 0) {
        print(".")
      }
      if (storeNodePHistory && iteration % (iterations / storeNodePHistIter).ceil == 0) {
        network.getNodes.foreach(node => node.probabilityHistory.add(node.posteriorProbability))
      }

      if (subnetworkSize.isEmpty)
        walkers.foreach(_.setGeneNumberVariable(4 + Random.nextInt(2))) // 4 Nodes => 3 interactions
      else
        walkers.foreach(_.setGeneNumberVariable(subnetworkSize.get))

      //walkers.foreach(_.setGeneNumberVariable(3 + StatUtils.getRandomPoisson(0.5))) // 0 - 60%, 1 - 30%, 2 - 10%)
      //walkers.foreach(_.setGeneNumberVariable(3))

      //.filterNot(networkManager.rankedGenes contains _)
      val callables = walkers.map(walker => {
        new Callable[Option[WalkerResult]] {
          override def call(): Option[WalkerResult] = {
        		  scoreWalker(walker)
          }
        }
      })

      val previous = System.nanoTime()

      val futures = callables.map(threadpool.submit(_))

      val after = System.nanoTime()

      //futures.map(x => println(x.get))
      val paths = futures.map(_.get.get)

      val genesScores = updateScores(paths)

      /*
         * Update converge iteration value if score converged
         */
      network.getNodes().foreach { node =>
        if (node.convergenceIteration == -1 && (node.posteriorProbability > 0.95 || node.posteriorProbability < 0.05)) {
          node.convergenceIteration = iteration
        }
      }

      if (iteration != iterations)
        evaporate(genesScores)

    }

    /*if (converged()) {
      println("Convergence reached [" + convergenceThreshold + "]")
    }*/

    if (!threadpool.isShutdown())
      threadpool.shutdown()
  }

}

