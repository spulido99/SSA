package be.cmpg.utils

import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.graph.Gene
import java.io.PrintWriter
import java.io.File

class MutualExclusivityPatternsManager(network: Network) {

  var allFoundSubnetworks = collection.mutable.Map[Network, Double]()
  var networksPerGeneMap = scala.collection.mutable.Map[String, Set[(Network, Double)]]()

  def addSubnetworks(subnetworks: scala.collection.Set[(scala.collection.Set[be.cmpg.graph.Interaction], Double)]) {

    subnetworks.foreach(subnetwork => {

      allFoundSubnetworks.put(new Network(subnetwork._1), subnetwork._2)

    })

    // This code was prohibitively slow
    /*    subnetworks.foreach(subnetwork => {
      allFoundSubnetworks.put(new Network(subnetwork._1), subnetwork._2)
    })*/

    //   val cutoff = subnetworks.size-(subnetworks.size/20)
    //   val sortedSubnetworks = subnetworks.toSeq.sortBy(_._2).toSet.drop(cutoff)
    //   // Only look at best 10% subnetworks
    //   
    //   sortedSubnetworks.map(element => element._1).flatten.foreach(interaction =>{
    //     patternNetworkInteractionsMap.put(interaction,patternNetworkInteractionsMap(interaction)+1)
    //   })

  }

  def calculateNbestSubnetworks {

    allFoundSubnetworks.filter(subnetwork => !subnetwork._2.isNaN).foreach(network => {

      val genesInNetwork = network._1.genes

      genesInNetwork.foreach(geneInNetwork => {
        if (networksPerGeneMap.get(geneInNetwork.name).isDefined) {

          networksPerGeneMap.put(geneInNetwork.name, networksPerGeneMap(geneInNetwork.name).+(network))
        } else {
          networksPerGeneMap.put(geneInNetwork.name, Set(network))
        }
      })

    })

    networksPerGeneMap = networksPerGeneMap.map(geneWithNetwork => {
      val gene = geneWithNetwork._1
      //     val filteredNetworkList = geneWithNetwork._2.toList.filter(network => geneWithNetwork._2.toList.filter(_._2 == network._2).head == network)
      val filteredNetworks = geneWithNetwork._2.toList.filter(network => geneWithNetwork._2.toList.filter(_._2 == network._2).head == network)
      val sortedNetworks = filteredNetworks.sortBy(_._2)
      // If a gene has more than 5 networks associated with it, only retain the 5 best.
      if (sortedNetworks.size > 5) {
        (gene, sortedNetworks.drop(sortedNetworks.size - 5).toSet)
      } else { (gene, sortedNetworks.toSet) }
    })
    allFoundSubnetworks = collection.mutable.Map[Network, Double]()
    
  }

  def printNetworkForSelectedGenes(selectedGenes: List[Gene], output: String, outputFolder: String) = {

    // Prints the interactions of the network which are between genes in a gene list.

    val patternNetworkInteractionsMap = collection.mutable.Map(network.interactions.map(interaction => (Interaction(interaction.from, interaction.to), 0)).toMap.toSeq: _*)

    val printM2 = new PrintWriter(new File(outputFolder + output + ".tsv"))

    patternNetworkInteractionsMap.filter(interaction => (interaction._1.listGenes.intersect(selectedGenes).size == 2)).foreach(interaction => {
      printM2.println(interaction._1.from.name + "\t" + interaction._1.to.name + "\t" + interaction._2)
    })

    printM2.close()
  }

  def returnNbestSubnetworks(selectedGenes: Set[Gene], output: String, outputFolder: String) {

    // Prints for every gene in the gene list the 5 best subnetworks in which that gene is involved.
    val printBestSubnetworks = new PrintWriter(new File(outputFolder + output + ".tsv"))

    networksPerGeneMap = networksPerGeneMap.filter(network => selectedGenes.contains(Gene(network._1)))
    
    networksPerGeneMap.foreach(geneWithNetwork => {
      val gene = geneWithNetwork._1
      val bestNetworks = geneWithNetwork._2
      bestNetworks.foreach(bestNetwork => { printBestSubnetworks.print(">" + gene + "," + bestNetwork._2 + "\n" + bestNetwork._1.toString() + "\n") })
    })

    printBestSubnetworks.close

    /*    val networksPerGeneMap = scala.collection.mutable.Map[String, Set[(Network, Double)]]()

    allFoundSubnetworks.filter(subnetwork => !subnetwork._2.isNaN).foreach(network => {

      val selectedGenesInNetwork = network._1.genes.intersect(selectedGenes)

      if (selectedGenesInNetwork.size > 0) {
        selectedGenesInNetwork.foreach(geneInNetwork => {
          if (networksPerGeneMap.get(geneInNetwork.name).isDefined) {

            networksPerGeneMap.put(geneInNetwork.name, networksPerGeneMap(geneInNetwork.name).+(network))
          } else {
            networksPerGeneMap.put(geneInNetwork.name, Set(network))
          }
        })

      }

    })

    networksPerGeneMap.foreach(geneWithNetwork => {
      val gene = geneWithNetwork._1
      val filteredNetworkList = geneWithNetwork._2.toList.filter(network => geneWithNetwork._2.toList.filter(_._2 == network._2).head == network)
      val sortedNetworks = geneWithNetwork._2.toList.sortBy(_._2)
      val NbestNetworksForGene = sortedNetworks.drop(sortedNetworks.size - N)
      NbestNetworksForGene.foreach(bestNetwork => printBestSubnetworks.print(">" + gene + "," + bestNetwork._2 + "\n" + bestNetwork._1.toString() + "\n"))
    })*/

    /*    val printBestSubnetworks = new PrintWriter(new File(output+".tsv"))
    selectedGenes.foreach(gene => {
      
    val subnetworksWithGene = allFoundSubnetworks.filter(network => network._1.genes.toList.intersect(selectedGenes).size > 0).toList
      
    val bestSubnetworksWithGene = subnetworksWithGene.sortBy(_._2).drop(subnetworksWithGene.size-5)
    bestSubnetworksWithGene.foreach(subnetwork => printBestSubnetworks.print(">"+gene+","+subnetwork._2+"\n"+subnetwork._1.toString()))
    })*/

  }

}