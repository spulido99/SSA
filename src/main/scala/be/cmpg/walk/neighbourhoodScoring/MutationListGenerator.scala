package be.cmpg.walk.neighbourhoodScoring
import be.cmpg.graph.Gene
import java.nio.file._
import java.io.File
import scala.io.Source

class MutationListGenerator(mutationfilesLocation: Path, referenceGenomeLocation: Path) {

  val mutationsDirectoryIterator = Files.newDirectoryStream(mutationfilesLocation).iterator()
  var allMutations = List[String]()

  def getExtendedMutationList: List[String] = {
    while (mutationsDirectoryIterator.hasNext()) {
      val file = mutationsDirectoryIterator.next().toString()
      val mutationsLines = Source.fromFile(new File(file)).getLines.toList
      val nameNumber = mutationsLines(0).split("\t").indexOf("locus_tag NC_000913K12MG1655 (CDS)")
      val synonymousnumber = mutationsLines(0).split("\t").indexOf("Non-synonymous")
      val locationInGenome = mutationsLines(0).split("\t").indexOf("Region")

      val mutationList = mutationsLines.drop(1).map(line => {
        val splitted = line.split("\t")
        if (splitted(nameNumber).contains("b")) {
          if (!(splitted(nameNumber).contains("\", \""))) {
            (splitted(nameNumber) + "\t" + splitted(synonymousnumber))
          } else {
            allMutations = allMutations.++(List(splitted(nameNumber).split("\", \"")(0) + "\t" + splitted(synonymousnumber)))
            (splitted(nameNumber).split("\", \"")(1) + "\t" + splitted(synonymousnumber))
          }
        } else {
          if (splitted(locationInGenome).contains("^")) {
            new IntergenicMutationAttributor(splitted(locationInGenome).split("\\^")(0).toInt, referenceGenomeLocation).Attribute + "\t" + splitted(synonymousnumber)
          } else if (splitted(locationInGenome).contains("..")) {
            new IntergenicMutationAttributor(splitted(locationInGenome).split("\\.\\.")(0).toInt, referenceGenomeLocation).Attribute + "\t" + splitted(synonymousnumber)
          } else {
            new IntergenicMutationAttributor(splitted(locationInGenome).toInt, referenceGenomeLocation).Attribute + "\t" + splitted(synonymousnumber)
          }
        }
      }).toList
      allMutations = allMutations.++(mutationList)
    }
    return allMutations.map(input => input.replaceAll("\"", ""))
  }
}