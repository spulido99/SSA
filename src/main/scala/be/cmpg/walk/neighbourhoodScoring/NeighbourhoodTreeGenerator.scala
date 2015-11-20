package be.cmpg.walk.neighbourhoodScoring
import be.cmpg.graph.Gene
import scala.util.control.Breaks._
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.expression.ExpressionNetworkManager

class NeighbourhoodTreeGenerator(startGene: Gene, depth: Int, networkManager: ExpressionNetworkManager) {

  private var visitedGenes = new scala.collection.mutable.HashSet[Gene]()
  private var leafs = new scala.collection.mutable.HashSet[Gene]()
  private var network = networkManager.getNetwork

  var score = networkManager.geneExpression(startGene);

  visitedGenes.add(startGene)
  leafs.add(startGene)

  def expand: (Double, scala.collection.mutable.HashSet[Gene]) = {
    var i = 0
    while (i < depth) {
      var tempLeafs = new scala.collection.mutable.HashSet[Gene]()

      leafs.foreach(leaf => {
        network.getAllInteractions(leaf).map(_.genes).flatten.foreach(newLeaf => {
          if (!visitedGenes.contains(newLeaf)) {
            score += networkManager.geneExpression(newLeaf);
            tempLeafs += newLeaf
            visitedGenes += newLeaf
          }
        })
      })
      leafs = tempLeafs
      i += 1
    }
    return (score / visitedGenes.size, visitedGenes)
  }
}