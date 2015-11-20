package be.cmpg.input

import org.specs2.mutable.Specification
import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import java.io.ByteArrayInputStream

@RunWith(classOf[JUnitRunner])
class TabHeaderFileSpecification extends Specification {

  "A TabHeaderFile" should {
    "be able to read a tab delimeted file with a header" in {
      val text = "id\tname\tfrom\na\tb\tc\nd\te\tf"

      val iterator = new TabHeaderFile(new ByteArrayInputStream(text.getBytes()))
      
      iterator.next.get("name").get must beEqualTo("b")
      iterator.next.get("name").get must beEqualTo("e")
        
    }
  }
}