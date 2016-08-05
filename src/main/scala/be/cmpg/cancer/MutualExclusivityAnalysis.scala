package be.cmpg.cancer

import scopt.OptionParser
import java.io.File
import be.cmpg.graph.Network
import java.util.HashSet
import be.cmpg.graph.Gene
import scala.collection.mutable.HashMap
import be.cmpg.walk.SubNetworkSelector
import scala.collection.JavaConversions._
import be.cmpg.utils.StatUtils
import java.io.FileWriter
import org.json.JSONObject
import org.json.JSONArray
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import be.cmpg.utils.MutualExclusivityPatternsManager
import scala.io.Source

object MutualExclusivityAnalysis extends App {

  val helper = new CancerHelper
  
  val parser = helper.getBasicArgParser("SSA.ME")
  
  // parser.parse returns Option[C]

  parser.parse(args, Config(

    /*
     * Default Values
     */
    iterations = 5000,
    reinforcement = 0.005,
    forgetfulness = 0.995
    
    /*,
    pvalIterations = 1000,
    pvalReinforcement = 0.005,
    pvalForgetfulness = 0.996
    */
    )) match {

    case Some(config) =>

      val (interactions, network, genePatientMatrix, all_samples, geneList) = helper.loadData(config,false)
      
      val mutualExclusivePatternsManager = new MutualExclusivityPatternsManager(network)
      
      println("Samples: " + all_samples.size)
      
      val networkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix,
        minimumSamplesAltered = config.minMutPerGene,
        pheromone = config.reinforcement,
        evaporation = config.forgetfulness,
        ranked = config.useRank,
        statistical = config.statistical
        )

      
      val walkers = helper.buildWalkers(geneList.view.toSet, networkManager)
      

      if (config.debug.isDefined)
        networkManager.debug = Some(config.debug.get.map { Gene(_) }.toSet, 20)
      //networkManager.debug = Some(Set(Gene("TP53")), 20)
      println("Subnetwork selectors: " + walkers.size)

      networkManager.run(config.iterations, geneList.view.toSet, config.processors, defwalkers = Some(walkers),mutualExclusivityPatternManager=Option(mutualExclusivePatternsManager),patterns=config.patterns)

      //val rankedGenes = networkManager.getRankedAllGenes().filter { g => networkManager.getPosteriorProbability(g) > 0.95 } .toList
      val rankedGenes = if (config.outputGenes == 0)
        networkManager.rankedGenes.toList
      else
        networkManager.getRankedAllGenes().take(config.outputGenes).toList

      val genesToReport = rankedGenes
      /*.map { g => network.getAllInteractions(g).map { i => i.genes }.flatten }.flatten
                          .filter { g => {
                            rankedGenes.contains(g) || network.getAllInteractions(g).map { i => i.genes }.flatten.toSet.intersect(rankedGenes).size > 2
                          }}*/

      /*val genesToReport = rankedGenes.toList
                          .map { g => network.getAllInteractions(g).map { i => i.genes }.flatten }.flatten
                          .groupBy(identity).mapValues(_.size)
                          .filter(e => rankedGenes.contains(e._1) || e._2 > 2)
                          .map(_._1)
                          .toSet*/

      println("\nGenes printed in network: " + rankedGenes.size + " " + rankedGenes.map(_.name))
      
      val malacardsGenes = if (config.otherGeneList==""){Set[Gene]()} else {Source.fromFile(config.otherGeneList).getLines.map(line => Gene(line.split("\t")(0))).toSet}
      
      MutualExclusivityPrintPattern.printPattern(config.outputPrefix, rankedGenes, networkManager, genePatientMatrix, config.outputFolder, config.inputFolder, otherPositiveGeneSetLists=malacardsGenes)

      if (config.patterns){
      println("\n calculating 5 best subnetworks per gene")
      
      mutualExclusivePatternsManager.returnNbestSubnetworks(rankedGenes.toSet,"5bestPatterns",config.outputFolder)
      }
      
      /**
       * P-value calculation
       */
      //val pvalueCalculator = new PValueCalculator(config.pvalIterations, genePatientMatrix, geneList)

    case None =>
      parser.showUsageAsError
    // arguments are bad, error message will have been displayed
  }

}