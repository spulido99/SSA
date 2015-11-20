package be.cmpg.walk.neighbourhoodScoring
import java.nio.file._
import scala.io.Source

class IntergenicMutationAttributor(MuationLocation: Int, referenceGenomeLocation: Path) {

  val referenceGenome = Source.fromFile(referenceGenomeLocation.toString()).getLines.toList
  val referenceGenes = referenceGenome.filter(entry => !(entry.startsWith("#")))
    .filter(entry => {
      val splitted = entry.split("\t")
      splitted(2).equals("gene")
    })
    .map(entry => {
      val splitted = entry.split("\t")
      (splitted(3) + "\t" + splitted(4), splitted(8).split("locus_tag=")(1).dropRight(splitted(8).split("locus_tag=")(1).length-5))
    }).toMap
    
  def Attribute:String = {
    var shortestDistance = 99999999
    var currentString = "None"
    
    referenceGenes.keys.toList.foreach(startpoint => {
      if (math.abs(startpoint.split("\t")(0).toInt - MuationLocation) < shortestDistance) { shortestDistance = math.abs(startpoint.split("\t")(0).toInt - MuationLocation)
        currentString = referenceGenes.get(startpoint).get
       }
    })
    
    referenceGenes.keys.toList.foreach(endpoint => {
      if (math.abs(endpoint.split("\t")(1).toInt - MuationLocation) < shortestDistance) { shortestDistance = math.abs(endpoint.split("\t")(1).toInt - MuationLocation) 
        currentString = referenceGenes.get(endpoint).get
       }
    })
    return currentString
  }

}