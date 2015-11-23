package be.cmpg.walk.neighbourhoodScoring

import java.nio.file.Paths
import be.cmpg.graph.Gene
import be.cmpg.graph.NetworkReader
import be.cmpg.graph.Network
import collection.immutable.ListMap
import java.io.FileWriter
import be.cmpg.graph.interaction.NodeCostNetworkManager
import scala.util.Random
import be.cmpg.walk.Path
import be.cmpg.utils.StatUtils
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.expression.ExpressionNetworkManager

object NeighbourhoodScanBySeedNode extends App {

  val writer = new FileWriter(raw"C:\Users\Bram\Documents\doctoraat\Projecten/PhAc/output.txt")
  val interactions = NetworkReader.fromFile(Paths.get("src/test/resources/be/cmpg/graph/network_small_connected.txt"))

  val selectedGenesSet = new FileWriter(raw"C:\Users\Bram\Documents\doctoraat\Projecten/PhAc/selected.txt")
  val goldenSet = new FileWriter(raw"C:\Users\Bram\Documents\doctoraat\Projecten/PhAc/goldenSet.txt")
  val possibleGenes = new FileWriter(raw"C:\Users\Bram\Documents\doctoraat\Projecten/PhAc/possibleGenes.txt")

  def simulateGeneScores(): (Set[Gene], ExpressionNetworkManager) = {
    /*
     * Select randomly the genes that are going to be important
     */

    var randomImportantGene = networkManager.getGenes.toList(Random.nextInt(networkManager.getGenes.size))
    //val randomImportantGene = Gene("STM4175")
    var resetcounter = 0
    var path = new Path(randomImportantGene)
    while (path.size < 7) {
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

    networkManager.getNetwork().getNodes().foreach( node => node.score = StatUtils.getRandomPoisson(if (importantGenes contains node.gene) importantGeneMean else nonImportantGeneMean))

    importantGenes.foreach(g => println(g + "\t" + networkManager.getNetwork().getNode(g).score))

    println("random scores generated")

    (importantGenes, networkManager)
  }
  
  val randomGeneratedInput = simulateGeneScores()
  val importantGenes = randomGeneratedInput._1
  println("importantGenes:" + importantGenes)
  val networkManager = randomGeneratedInput._2

  importantGenes.foreach(gene => goldenSet.write(gene.name + "\n"))
  goldenSet.close

  networkManager.getNetwork().getNodes().foreach { node => possibleGenes.write(node.gene.name + "\n") }
  possibleGenes.close

  val mutationsScore = networkManager.getNetwork().getNodes().map(node => (node.gene.name, node.score.toString)).toMap
  var relevanceMap = new scala.collection.mutable.HashMap[String, Double]()

  println("mutationScores:" + mutationsScore)

  networkManager.getGenes.foreach(mutatedGene => {
    val MutationDisRelevance = new NeighbourhoodTreeGenerator(mutatedGene, 2, networkManager).expand
    relevanceMap.put(mutatedGene.name, MutationDisRelevance._1)
  })

  val orderedMap = ListMap(relevanceMap.toList.sortBy { _._2 }: _*)

  writer.write("gene_name" + "\t" + "Simulated Neighbourhood_score" + "\t" + "individual simulated score" + "\t" + "Number of non-synonymous mutations" + "\t" + "number of synonymous mutations" + "\n")

  orderedMap.foreach(entry => {
    writer.write(entry._1 + "\t" + entry._2 + "\t" + mutationsScore.get(entry._1).get + "\n")
  })
  writer.close
  println("orderedResults" + orderedMap)
  val selectedGenesWithValues = orderedMap.map(input => input._1 + "\t" + input._2.toString).toList.dropRight(orderedMap.size - 7)

  val selectedGenes = selectedGenesWithValues.map(element => element.split("\t")(0)).toSet

  selectedGenes.foreach(gene => selectedGenesSet.write(gene))

  selectedGenesSet.close
}