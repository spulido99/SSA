package be.cmpg.cancer.simulations

import be.cmpg.cancer.CancerHelper
import scala.util.Random
import be.cmpg.walk.Path
import be.cmpg.cancer.MutualExclusivityNetworkManager
import be.cmpg.graph.Network
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.JavaConversions._
import java.util.HashMap
import be.cmpg.cancer.Polimorphism
import be.cmpg.graph.Gene
import java.util.HashSet
import be.cmpg.simulatedDataAnalysis.ROCPlotValuesCreator
import be.cmpg.graph.Interaction
import java.io.FileWriter
import be.cmpg.cancer.MutualExclusivityPrintPattern
import be.cmpg.cancer.PolimorphismKey
import be.cmpg.cancer.PolimorphismKey
import be.cmpg.cancer.Polimorphism
import java.io.PrintWriter

object MutualExclusivitySimulations extends App {

  /* Helpers */
  val helper = new CancerHelper()
  val random = new Random(System.nanoTime())
  val development = args.length == 0

  /* simulated pattern parameters */
  val meSamplesPercentaje = 0.30
  val subnetworkSize = 20

  /* simulated pattern parameters for Convergence Testing */
  //println( ">>>> CONVERGENCE PARAMETERS USED <<<<")
  //val meSamplesPercentaje = 0.80
  //val subnetworkSize = 20

  /* method parameters */

  val noiseLevels = if (development) List(0.0) else List(0.0, 0.1, 0.25, 0.50)

  /* simulation parameters */
  val repeats = if (development) 1 else args(0).toInt

  val subnetworkSizePerformance = true
  val runTimeEstimation = false
  val runConvergence = false
  val runParametersEstimation = false
  val runRobustisity = false

  val reinforcement = 0.005
  val forgetfulness = 0.995

  val seedThreshold = 10
  val iterations = 1000

  println("*************************************************")
  if (development) {
    println("*********** DEVELOPMENT *************************")
  } else {
    println("*************** SERVER **************************")
  }
  println("*************************************************")

  /* Constants */
  val networkName = "HT_hiII14"
  val interactions = helper.loadReferenceNetwork(networkName)
  val baseNetwork = new Network(interactions)
  println("Base Network size: " + baseNetwork.interactions.size)
  val genes = baseNetwork.genes.toList
  val mutQtyBySample = new CSVReader(new FileReader("MutualExclusivityStatistics.txt"), '\t')
    .readAll()
    .map { f => (f(0), f(1).toInt) }

  /*
   * Subnetwork size performance
   */
  
  if (subnetworkSizePerformance) {
    val writerSubnetworksize = new FileWriter(networkName + "simulation_output.subnetworksize.txt")

    writerSubnetworksize.write("NoiseLevel\tRep\tN\tTP\tFP\tTN\tFN\n")

    for (r <- 0 until 10) {

      val causalNetwork = getRandomCausalNetwork(baseNetwork)
      /*val causalNetwork = new Network(Set(
        Interaction(Gene("CCDC60"), Gene("MTUS2")), 
        Interaction(Gene("MTUS2"), Gene("SNHG11")), 
        Interaction(Gene("LZTS2"), Gene("SNHG11")), 
        Interaction(Gene("CCNG1"), Gene("LZTS2"))))
      */

      // (gene, sample) -> Mutation
      /*
       * Test finding ME pattern
       */
      val genePatientMatrix: Map[PolimorphismKey, Polimorphism] = buildRandomGenePatientMatrix ++ buildMEGenePatientMatrix(causalNetwork)
      /*
       * Test just random data
       */
      //val genePatientMatrix = buildRandomGenePatientMatrix

      val geneMutCounts = genePatientMatrix.groupBy(e => e._1.gene).map(e => (e._1, e._2.size))

      for (size <- 15.to(30, 5)) {//1 to 10) {
        println("*************************************************")
        println(">>>> Subnetwork Size: " + size + " repeat: " + r)
        println("*************************************************")

        println("ME genes: ")
        causalNetwork.genes.foreach { g => println(g.name + " > " + geneMutCounts.getOrElse(g, 0)) }

        val copyNetwork = new Network(baseNetwork.interactions)

        val geneList = geneMutCounts.filter(e => e._2 >= seedThreshold).map(_._1).toSet

        println("Subnetwork Sizes simulation")
        runSimulation(genePatientMatrix, causalNetwork, size, r, writerSubnetworksize, copyNetwork, geneList, Some(size))

      }
    }

    writerSubnetworksize.close()
  }
  
  /*
   * Time stimation
   */

  if (runTimeEstimation) {

    val writerRemove = new PrintWriter(networkName + "simulation_output.processingtime.txt")
    writerRemove.println("Rep\tSeedQty\tTimeNano")

    val threadpool = java.util.concurrent.Executors.newFixedThreadPool(10)

    for (r <- 0 until 10) {

      val causalNetwork = getRandomCausalNetwork(baseNetwork)
      /*val causalNetwork = new Network(Set(
        Interaction(Gene("CCDC60"), Gene("MTUS2")), 
        Interaction(Gene("MTUS2"), Gene("SNHG11")), 
        Interaction(Gene("LZTS2"), Gene("SNHG11")), 
        Interaction(Gene("CCNG1"), Gene("LZTS2"))))
      */

      // (gene, sample) -> Mutation
      /*
       * Test finding ME pattern
       */
      val genePatientMatrix: Map[PolimorphismKey, Polimorphism] = buildRandomGenePatientMatrix ++ buildMEGenePatientMatrix(causalNetwork)

      val geneMutCounts = genePatientMatrix.groupBy(e => e._1.gene).map(e => (e._1, e._2.size)).toList.sortBy(-_._2)

      for (percent <- 0.0 until (1.0, 0.1)) {
        val numberOfSeedGenes = (geneMutCounts.size * percent).toInt max 1

        threadpool.submit(new Runnable() {
          def run() {

            val copyNetwork = new Network(baseNetwork.interactions)

            val networkManager = new MutualExclusivityNetworkManager(
              network = copyNetwork,
              genePatientMatrix = genePatientMatrix,
              minimumSamplesAltered = 3,
              pheromone = 0.005,
              evaporation = 0.99,
              ranked = true,
              minProb = 0.0)

            val geneList = geneMutCounts.take(numberOfSeedGenes).map(_._1).toSet
            println(r + "\t" + numberOfSeedGenes)
            val start = System.nanoTime()
            networkManager.run(5000, geneList, 1, storeNodePHistory = false)
            val end = System.nanoTime()
            println(r + "\t" + numberOfSeedGenes + "\t" + (end - start))
            writerRemove.println(r + "\t" + numberOfSeedGenes + "\t" + (end - start))
          }
        })

      }

    }
    threadpool.shutdown
    writerRemove.close

  }

  /*
   * Convergence
   */

  if (runConvergence) {

    val causalNetwork = getRandomCausalNetwork(baseNetwork)
    val genePatientMatrix: Map[PolimorphismKey, Polimorphism] = buildRandomGenePatientMatrix ++ buildMEGenePatientMatrix(causalNetwork)

    val geneMutCounts = genePatientMatrix.groupBy(e => e._1.gene).map(e => (e._1, e._2.size))

    val geneList = geneMutCounts.filter(e => e._2 >= seedThreshold).map(_._1).toSet

    val copyNetwork = new Network(baseNetwork.interactions)

    val networkManager = new MutualExclusivityNetworkManager(
      network = copyNetwork,
      genePatientMatrix = genePatientMatrix,
      minimumSamplesAltered = 3,
      pheromone = 0.05,
      evaporation = 0.990,
      ranked = true,
      minProb = 0.0)

    //networkManager.debugForGenes(causalNetwork.genes, 20)
    networkManager.run(50000, geneList, Runtime.getRuntime().availableProcessors() / 2, storeNodePHistory = true, storeNodePHistIter = 500)

    val historicPWriter = new PrintWriter("genesHistoricPs.out")
    networkManager.getRankedAllGenes().foreach(g => {
      historicPWriter.print(g.name)
      val node = networkManager.getNetwork.getNode(g)
      node.probabilityHistory.foreach(p => historicPWriter.print("\t" + (p * 100).toInt))
      historicPWriter.println()
    })
    historicPWriter.close()
  }

  /*
   * Robustisity
   */

  if (runRobustisity) {

    val writerRemove = new FileWriter(networkName + "simulation_output.remove.txt")
    val writerAdded = new FileWriter(networkName + "simulation_output.added.txt")

    writerRemove.write("NoiseLevel\tRep\tN\tTP\tFP\tTN\tFN\n")
    writerAdded.write("NoiseLevel\tRep\tN\tTP\tFP\tTN\tFN\n")

    for (r <- 0 until repeats) {

      val causalNetwork = getRandomCausalNetwork(baseNetwork)
      /*val causalNetwork = new Network(Set(
        Interaction(Gene("CCDC60"), Gene("MTUS2")), 
        Interaction(Gene("MTUS2"), Gene("SNHG11")), 
        Interaction(Gene("LZTS2"), Gene("SNHG11")), 
        Interaction(Gene("CCNG1"), Gene("LZTS2"))))
      */

      // (gene, sample) -> Mutation
      /*
       * Test finding ME pattern
       */
      val genePatientMatrix: Map[PolimorphismKey, Polimorphism] = buildRandomGenePatientMatrix ++ buildMEGenePatientMatrix(causalNetwork)
      /*
       * Test just random data
       */
      //val genePatientMatrix = buildRandomGenePatientMatrix

      val geneMutCounts = genePatientMatrix.groupBy(e => e._1.gene).map(e => (e._1, e._2.size))

      for (noise <- noiseLevels) {
        println("*************************************************")
        println(">>>> Data Noise: " + noise + " repeat: " + r)
        println("*************************************************")

        println("ME genes: ")
        causalNetwork.genes.foreach { g => println(g.name + " > " + geneMutCounts.getOrElse(g, 0)) }

        val noisyNetwork_removed = new Network(baseNetwork.interactions.filter(i => random.nextDouble() >= noise))

        val interactionsToAdd = (0L until (baseNetwork.interactions.size * noise).round).map { i =>
          Interaction(genes(random.nextInt(genes.size)), genes(random.nextInt(genes.size)))
        }
        val noisyNetwork_added = new Network(baseNetwork.interactions ++ interactionsToAdd)

        val geneList_removed = geneMutCounts.filter(e => e._2 >= seedThreshold && noisyNetwork_removed.genes.contains(e._1)).map(_._1).toSet

        val geneList_added = geneMutCounts.filter(e => e._2 >= seedThreshold).map(_._1).toSet

        println("Removed edges simulation")
        runSimulation(genePatientMatrix, causalNetwork, noise, r, writerRemove, noisyNetwork_removed, geneList_removed)
        println("Added edges simulation")
        runSimulation(genePatientMatrix, causalNetwork, noise, r, writerAdded, noisyNetwork_added, geneList_added)

      }
    }

    writerRemove.close()
    writerAdded.close()

  }
  /*
   * Parameter sensitivity
   */

  if (runParametersEstimation) {

    val rs = (1 to 10).map(0.0010 * _)
    val fs = rs.map(1.0 - _)

    println("*****************")
    println(">> Parameters <<")
    println(">> Reinforncements: " + rs)
    println(">> Forgetfulness: " + fs)

    val writerParameters = new FileWriter(networkName + "simulation_output.parameters.txt")

    writerParameters.write("Rep\tReinforcement\tForgetfullness\tN\tTP\tFP\tTN\tFN\n")

    val reps = if (development) 1 to 1 else 1 to 10

    for (rep <- reps) {

      val causalNetwork = getRandomCausalNetwork(baseNetwork)
      val genePatientMatrix: Map[PolimorphismKey, Polimorphism] = buildRandomGenePatientMatrix ++ buildMEGenePatientMatrix(causalNetwork)

      val geneMutCounts = genePatientMatrix.groupBy(e => e._1.gene).map(e => (e._1, e._2.size))

      //for (p <- parameters) {
      for (r <- rs) {
        for (f <- fs) {

          val p = Param(r, f)

          println("*******************************")
          println(">>> Repetition: " + rep + " Reinforcemnt: " + p.reinforcement + " Forguetfulness: " + p.forgetfulness)
          println("*******************************")

          val geneList = geneMutCounts.filter(e => e._2 >= seedThreshold).map(_._1).toSet

          println("Number of Walkers: " + geneList.size)
          println(causalNetwork.interactions)

          val copyNetwork = new Network(baseNetwork.interactions)

          val networkManager = new MutualExclusivityNetworkManager(
            network = copyNetwork,
            genePatientMatrix = genePatientMatrix,
            minimumSamplesAltered = 3,
            pheromone = p.reinforcement,
            evaporation = p.forgetfulness,
            ranked = true)

          //networkManager.debugForGenes(causalNetwork.genes, 20)
          networkManager.run(iterations, geneList, Runtime.getRuntime().availableProcessors() / 2)

          val rankedGenes = networkManager.getRankedAllGenes().filter { g => networkManager.getPosteriorProbability(g) > 0.95 }.toList
          println("Genes selected: " + rankedGenes.size)
          println(rankedGenes)

          if (development) {
            MutualExclusivityPrintPattern.printPattern("sim", rankedGenes, networkManager, genePatientMatrix)
          }

          val rankedGeneByPhac = networkManager.getRankedAllGenes().map(_.name).toList;

          val importantGenes = causalNetwork.genes.map { g => g.name }.view.toSet
          val allGenes = baseNetwork.genes.map { g => g.name }.view.toSet

          val phacSelected = new HashSet[String]
          for (n <- 0 to baseNetwork.genes.size - 1) {
            val gene = rankedGeneByPhac(n)
            phacSelected += gene
            val phacROCValues = new ROCPlotValuesCreator(phacSelected.view.toSet, importantGenes, allGenes).get
            if (importantGenes contains gene)
              println(gene + "\t" + n + "\t" + phacROCValues._TP + "\t" + phacROCValues._FP + "\t" + phacROCValues._TN + "\t" + phacROCValues._FN)

            writerParameters.write(rep + "\t" + p.reinforcement + "\t" + p.forgetfulness + "\t" + n + "\t" + phacROCValues._TP + "\t" + phacROCValues._FP + "\t" + phacROCValues._TN + "\t" + phacROCValues._FN + "\n")
          }
        }
      }
    }
    writerParameters.close()
  }

  /*
   ************************************
   * Methods
   ************************************
   */

  def runSimulation(genePatientMatrix: Map[PolimorphismKey, Polimorphism], causalNetwork: Network, noise: Double, r: Int, writer: FileWriter, noisyNetwork: Network, geneList: Set[Gene], subnetworkSize:Option[Int]=None) = {
    println("Number of Walkers: " + geneList.size)
    println("Causal Network: " + causalNetwork.interactions)
    println("Noisy Network Size: " + noisyNetwork.interactions.size)

    val networkManager = new MutualExclusivityNetworkManager(
      network = noisyNetwork,
      genePatientMatrix = genePatientMatrix,
      minimumSamplesAltered = 3,
      pheromone = reinforcement,
      evaporation = forgetfulness,
      ranked = true)

    //networkManager.debugForGenes(causalNetwork.genes, 20)
    networkManager.run(iterations, geneList, Runtime.getRuntime().availableProcessors() / 2)

    val rankedGenes = networkManager.getRankedAllGenes().filter { g => networkManager.getPosteriorProbability(g) > 0.95 }.toList
    println("Genes selected: " + rankedGenes.size)
    println(rankedGenes)

    if (development) {
      MutualExclusivityPrintPattern.printPattern("sim", rankedGenes, networkManager, genePatientMatrix)
    }

    val rankedGeneByPhac = networkManager.getRankedAllGenes().map(_.name).toList;

    val importantGenes = causalNetwork.genes.map { g => g.name }.view.toSet
    val allGenes = noisyNetwork.genes.map { g => g.name }.view.toSet

    val phacSelected = new HashSet[String]
    for (n <- 0 to noisyNetwork.genes.size - 1) {
      val gene = rankedGeneByPhac(n)
      phacSelected += gene
      val phacROCValues = new ROCPlotValuesCreator(phacSelected.view.toSet, importantGenes, allGenes).get
      if (importantGenes contains gene)
        println(gene + "\t" + n + "\t" + phacROCValues._TP + "\t" + phacROCValues._FP + "\t" + phacROCValues._TN + "\t" + phacROCValues._FN)

      writer.write(noise + "\t" + r + "\t" + n + "\t" + phacROCValues._TP + "\t" + phacROCValues._FP + "\t" + phacROCValues._TN + "\t" + phacROCValues._FN + "\n")
    }
  }

  def getRandomCausalNetwork(network: Network) = {

    var randomImportantGene = genes(random.nextInt(genes.size))
    //val randomImportantGene = Gene("STM4175")
    var resetcounter = 0
    var path = new Path(randomImportantGene)
    while (path.size < subnetworkSize) {
      if (random.nextDouble <= 1 / subnetworkSize) {
        path.restart
      } else {
        val possibleInteractions = network.getOutgoingInteractions(path.currentEndpoint).filter(interaction => path.canTakeInteraction(interaction)).toList
        if (!possibleInteractions.isEmpty) {
          val next = possibleInteractions(random.nextInt(possibleInteractions.size))
          path.expand(next)
        } else {
          resetcounter = resetcounter + 1
          if (resetcounter > 20) {
            randomImportantGene = genes(random.nextInt(genes.size))
            path = new Path(randomImportantGene)
            resetcounter = 0
          }
          path.reset
        }
      }
    }

    new Network(path.getVisitedInteractions().toSet)
  }

  def buildMEGenePatientMatrix(causalNetwork: Network) = {

    // randomly select 30% of samples to be the ones with the mutually exclusive pattern
    val meSamples = mutQtyBySample.filter(s => random.nextDouble() < meSamplesPercentaje)
    /*List(
    ("TCGA-A2-A0SW",  62),
    ("TCGA-E2-A15O",  46),
    ("TCGA-A2-A1G6",  7))*/

    println("Samples in ME: " + meSamples.size)

    val samplesToAssign = meSamples.toBuffer

    val meGenes = random.shuffle(causalNetwork.genes.toList)

    // Some Genes that are not mutated but are part of the pattern
    //val noMutableMEGenes = meGenes.filter { g => random.nextDouble() < 0.2 }

    val meGenePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]

    var index = 0
    while (!samplesToAssign.isEmpty) {
      //if (random.nextDouble() < 1.0 / (index + 1)){ //math.min((index + 2.0), 15.0)) {
      //if (random.nextDouble() < 1.0 - (index/subnetworkSize.doubleValue)){
      if (random.nextDouble() < 1.0 / math.min((index + 1.0) * 3, 30.0)) {
        // assign mutation in ME gene
        val gene = meGenes(index)
        if (true) { //}!noMutableMEGenes.contains(gene)) {
          val sample = samplesToAssign(0)._1
          meGenePatientMatrix.put(PolimorphismKey(gene, sample), Polimorphism(gene.name))

          // have random mutation in the other genes?
          /*
           * 20% prob that the network is mutated somewhere else.
           * Number as a random because its probably to be mutated more than once
           */
          meGenes.filter { g => random.nextDouble() < 0.2 / subnetworkSize }.foreach { g => meGenePatientMatrix.put(PolimorphismKey(g, sample), Polimorphism(g.name)) }

          samplesToAssign.remove(0)
        }
      }
      index = (index + 1) % meGenes.size
    }

    meGenePatientMatrix
  }

  def buildRandomGenePatientMatrix = mutQtyBySample.map { m =>
    {
      (0 until m._2).map { i =>
        {
          val randomGene = genes(random.nextInt(genes.size))
          val key = PolimorphismKey(randomGene, m._1)
          val value = Polimorphism(randomGene.name)
          (key, value)
        }
      }
    }
  }.flatten.toMap

}
case class Param(reinforcement: Double, forgetfulness: Double)