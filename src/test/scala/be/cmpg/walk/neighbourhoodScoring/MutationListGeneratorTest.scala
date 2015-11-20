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
object MutationListGeneratortest extends Specification {

  "a MutationListGenerator" should {
    "be able to generate a list of mutations and its synonymity state based on a number of mutation files in a directory" in {
      val mutationDirectoryPath = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/mutations")
      
      val geneReferenceFilePath = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/referenceGenomes/EcoliK12MG1655.gff")
      
      val mutationList = new MutationListGenerator(mutationDirectoryPath,geneReferenceFilePath).getExtendedMutationList
      
      mutationList must be equalTo List("b001"+"\t"+"Yes","b009"+"\t"+"Yes","b001"+"\t"+"No","b001"+"\t"+"-","b006"+"\t"+"Yes","b003"+"\t"+"Yes","b003"+"\t"+"No","b003"+"\t"+"-")
    }
    
  }

}