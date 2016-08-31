package be.cmpg.walk.neighbourhoodScoring

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import java.nio.file.Paths
import be.cmpg.graph.Gene
import be.cmpg.graph.NetworkReader
import be.cmpg.graph.Network
import collection.immutable.ListMap
import java.io.FileWriter
import be.cmpg.expression.ExpressionNetworkManager
import be.cmpg.utils.weightByFlatInitialProbability

@RunWith(classOf[JUnitRunner])
class UniformNetworkTest extends Specification {

  "A NeighbourhoodTreeGenerator" should {
    "give back a acore of 5 if every gene in the network has a score of 5" in {

      val mutationDirectoryPath = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/uniformScoresMutations/")

      val geneReferenceFilePath = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/referenceGenomes/EcoliK12MG1655.gff")

      val mutationList = new MutationListGenerator(mutationDirectoryPath, geneReferenceFilePath).getExtendedMutationList

      val mutationScoringMap = new MutationScoresGenerator(mutationList).getScoresCount.map(mutScore => (Gene(mutScore._1),mutScore._2.split("\t")(0).toDouble)).toMap

      val network = new Network(NetworkReader.fromFile(raw"src\test\resources\be\cmpg\NeighbourhoodScoring\graph\Network.txt"))

      val startGene = new Gene("b001")

      val searchdepth = 3

      val networkManager = new ExpressionNetworkManager(network, weightingScheme = new weightByFlatInitialProbability(network, 0.5))

      networkManager.geneExpression = mutationScoringMap
      
      val normalizedScore = new NeighbourhoodTreeGenerator(startGene, searchdepth, networkManager).expand

      normalizedScore._1 must be equalTo 5

    }
  }
}