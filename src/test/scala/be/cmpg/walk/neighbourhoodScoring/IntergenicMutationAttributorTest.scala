package be.cmpg.walk.neighbourhoodScoring
import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import java.nio.file._

@RunWith(classOf[JUnitRunner])
class IntergenicMutationAttributorTest extends Specification {

  "an IntergenicMutationAttributorTest" should{
    "be able to assign the closest gene to an intergenic mutation" in{
      
      val referenceGenome = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/referenceGenomes/EcoliK12MG1655.gff")
      
      val mutationsDirectory = Paths.get("src/test/resources/be/cmpg/NeighbourhoodScoring/IntergenicMutations")
      // Note that also deletion, insertions and mutations in multiple genes at the same time are tested here
        new MutationListGenerator(mutationsDirectory,referenceGenome).getExtendedMutationList must be equalTo(List("b0556"+"\t"+"Yes","b0005"+"\t"+"Yes","b0023" +"\t"+"Yes","b4510"+"\t"+"Yes","b0023"+"\t"+"Yes"))
    }
  }
  
}