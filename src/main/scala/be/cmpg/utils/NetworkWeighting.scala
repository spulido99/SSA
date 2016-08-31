package be.cmpg.utils
import be.cmpg.graph.Network
import org.apache.commons.math3.distribution.ExponentialDistribution

abstract class NetworkWeighting {

  def initialize: Network
}

// Weigh initial probability of genes based on the exponential distribution of the degrees. Genes with high degree are downweighted
class weightByDegreeExponentialdistribution(network: Network, initialProb: Double) extends NetworkWeighting {

  def initialize(): Network = {
    val meanDegreeOfNetwork = network.genes.toList.map(gene => { network.getInDegree(gene) + network.getOutDegree(gene) }).sum.toDouble / (network.genes.toList.map(gene => { network.getInDegree(gene) + network.getOutDegree(gene) }).size.toDouble)

    val exponentialDistribution = new ExponentialDistribution(meanDegreeOfNetwork)

    network.getNodes().foreach { node =>
      val probability = 1 - exponentialDistribution.cumulativeProbability(network.getInDegree(node.gene) + network.getOutDegree(node.gene))
      node.posteriorProbability = initialProb * probability
    }

    network

  }
}

// Assign a flat initial probability to all genes in the network
class weightByFlatInitialProbability(network: Network, initialProb: Double) extends NetworkWeighting {

  def initialize(): Network = {
    network.getNodes().foreach { node =>
      node.posteriorProbability = initialProb
    }

    network
  }
}

// Assign a specific probability to a list of genes in the network. The other genes get the initialprobability. This method is mainly used for testing.
class testWeighting(network: Network, initialProb: Double,nodeNames:List[String],specificProb:Double) extends NetworkWeighting {
  def initialize(): Network = {
    network.getNodes().foreach { node =>
      if (nodeNames.contains(node.gene.name)) { node.posteriorProbability = specificProb }
      else{
      node.posteriorProbability = initialProb}
    }

    network
  }
}