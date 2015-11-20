package be.cmpg.input

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import scala.io.Source
import java.io.ByteArrayInputStream

@RunWith(classOf[JUnitRunner])
class ScoreFileSpecification extends Specification {

  "A ScoreFile" should {

    val content = """
biofilm_BF_14h_vs_BF_24h,CHIX,0.6491000000000007
biofilm_BF_7h_vs_BF_10h,CSRB,0.8239999999999998
biofilm_BF_10h_vs_BF_14h,CSRB,0.036999999999999034
biofilm_BF_14h_vs_BF_24h,CSRB,0.04600000000000115
biofilm_BF_7h_vs_BF_10h,CSRC,1.5030000000000001
biofilm_BF_10h_vs_BF_14h,CSRC,0.2029999999999994"""

    val scoreFile = new ScoreFile(new ByteArrayInputStream(content.getBytes()))

    "be able to read the different conditions in the file" in {
      scoreFile.getConditions must contain( allOf(Set("biofilm_BF_14h_vs_BF_24h", "biofilm_BF_10h_vs_BF_14h", "biofilm_BF_7h_vs_BF_10h")))
    }
    
    "be able to parse all the values" in {
      scoreFile.getValues must contain(("biofilm_BF_14h_vs_BF_24h","CHIX",0.6491000000000007))
      scoreFile.getValues must have size(6)
    }
    
    "be able to get all score pairs for a single condition" in{
      scoreFile.getScoresForCondition("biofilm_BF_14h_vs_BF_24h") must contain (allOf(("CSRB",0.04600000000000115), ("CHIX",0.6491000000000007)))
    }
  }
}