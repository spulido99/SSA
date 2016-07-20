package be.cmpg.utils

import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.graph.Gene
import java.io.PrintWriter
import java.io.File

class MutualExclusivityPatternsManager(network: Network) {

  val allFoundSubnetworks = collection.mutable.Map[Network, Double]()

  def addSubnetworks(subnetworks: scala.collection.Set[(scala.collection.Set[be.cmpg.graph.Interaction], Double)]) {

    subnetworks.foreach(subnetwork => {
      allFoundSubnetworks.put(new Network(subnetwork._1), subnetwork._2)
    })

    //   val cutoff = subnetworks.size-(subnetworks.size/20)
    //   val sortedSubnetworks = subnetworks.toSeq.sortBy(_._2).toSet.drop(cutoff)
    //   // Only look at best 10% subnetworks
    //   
    //   sortedSubnetworks.map(element => element._1).flatten.foreach(interaction =>{
    //     patternNetworkInteractionsMap.put(interaction,patternNetworkInteractionsMap(interaction)+1)
    //   })

  }

  def printNetworkForSelectedGenes(selectedGenes: List[Gene], output: String) = {

    // Prints the interactions of the network which are between genes in a gene list.

    val patternNetworkInteractionsMap = collection.mutable.Map(network.interactions.map(interaction => (Interaction(interaction.from, interaction.to), 0)).toMap.toSeq: _*)

    val printM2 = new PrintWriter(new File(output + ".tsv"))

    patternNetworkInteractionsMap.filter(interaction => (interaction._1.listGenes.intersect(selectedGenes).size == 2)).foreach(interaction => {
      printM2.println(interaction._1.from.name + "\t" + interaction._1.to.name + "\t" + interaction._2)
    })

    printM2.close()
  }

  def returnNbestSubnetworks(selectedGenes: Set[Gene], N: Int, output: String) {

    // Prints for every gene in the gene list the 5 best subnetworks in which that gene is involved.
    val printBestSubnetworks = new PrintWriter(new File(output+".tsv"))
    
    val networksPerGeneMap = scala.collection.mutable.Map[String, Set[(Network, Double)]]()
    
    allFoundSubnetworks.filter(subnetwork => !subnetwork._2.isNaN).foreach(network => {

      val selectedGenesInNetwork = network._1.genes.intersect(selectedGenes)

      if (selectedGenesInNetwork.size > 0) {
        selectedGenesInNetwork.foreach(geneInNetwork => {
          if (networksPerGeneMap.get(geneInNetwork.name).isDefined){
            networksPerGeneMap.put(geneInNetwork.name,networksPerGeneMap(geneInNetwork.name).+(network))
          }
          else{
          networksPerGeneMap.put(geneInNetwork.name, Set(network))}
        })

      }
      
    })
    
    val bestNetworksMap = networksPerGeneMap.foreach(geneWithNetwork =>{
      val gene = geneWithNetwork._1
      val sortedNetworks = geneWithNetwork._2.toList.sortBy(_._2)
      val NbestNetworksForGene = sortedNetworks.drop(sortedNetworks.size-N)
      NbestNetworksForGene.foreach(bestNetwork => printBestSubnetworks.print(">"+gene+","+bestNetwork._2+"\n"+bestNetwork._1.toString()+"\n"))
    })
    
    printBestSubnetworks.close

    /*    val printBestSubnetworks = new PrintWriter(new File(output+".tsv"))
    selectedGenes.foreach(gene => {
      
    val subnetworksWithGene = allFoundSubnetworks.filter(network => network._1.genes.toList.intersect(selectedGenes).size > 0).toList
      
    val bestSubnetworksWithGene = subnetworksWithGene.sortBy(_._2).drop(subnetworksWithGene.size-5)
    bestSubnetworksWithGene.foreach(subnetwork => printBestSubnetworks.print(">"+gene+","+subnetwork._2+"\n"+subnetwork._1.toString()))
    })*/

  }

}