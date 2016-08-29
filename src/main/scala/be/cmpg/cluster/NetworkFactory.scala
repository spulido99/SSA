package be.cmpg.cluster

import java.io.File
import be.cmpg.graph.Interaction
import scala.io.Source
import be.cmpg.graph.Gene

object NetworkFactory {
  def loadNetwork(networkFile: File): Set[Interaction] = {
    Source.fromFile(networkFile).getLines.map(line => {
      val splitted = line.split(",")
      Interaction(new Gene(splitted(0)), new Gene(splitted(1)))
    }).toSet
  }
}