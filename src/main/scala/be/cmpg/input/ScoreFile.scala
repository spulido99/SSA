package be.cmpg.input

import java.io.InputStream
import scala.io.Source

class ScoreFile(stream: InputStream) {

  private def getLines =
    Source.fromInputStream(stream)
      .getLines
      .filter(!_.isEmpty())
      .filter(!_.startsWith("#"))

  lazy val getConditions: Set[String] = {
    getValues.map(_._1).toSet
  }

  lazy val getValues: List[(String, String, Double)] = {
    getLines.map((line) => {
      val splitted = line.split(",")
      (splitted(0), splitted(1), splitted(2).toDouble)
    }).toList
  }
  
  def getScoresForCondition(condition: String): List[(String, Double)] = {
    getValues.filter(_._1 == condition).map(x => (x._2, x._3))
  }
}