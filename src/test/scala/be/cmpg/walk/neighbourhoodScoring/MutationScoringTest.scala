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

@RunWith(classOf[JUnitRunner])
class MutationScoringTest extends Specification {

  "a MutationScoringGenerator" should {
    "generate a map of genes with their associated scores" in {

      val mutationDirectoryPath = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/mutations/")

      val geneReferenceFilePath = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/referenceGenomes/EcoliK12MG1655.gff")

      val mutationList = new MutationListGenerator(mutationDirectoryPath, geneReferenceFilePath).getExtendedMutationList

      val mutationScoringMap = new MutationScoresGenerator(mutationList).getScoresQuotient

      mutationScoringMap must be equalTo (scala.collection.mutable.HashMap[String, String](("b009", "1.0\t1\t0"), ("b003", (2d / 3d).toString+"\t2\t1"), ("b006", "1.0\t1\t0"), ("b001", (2d / 3d).toString+"\t2\t1")))
    }
  }
}