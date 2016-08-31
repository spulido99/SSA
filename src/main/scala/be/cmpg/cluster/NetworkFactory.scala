package be.cmpg.cluster

import java.io.File
import be.cmpg.graph.Interaction
import scala.io.Source
import be.cmpg.graph.Gene
import be.cmpg.graph.Network

object NetworkFactory {
  def loadNetwork(networkFile: File): Network = {
    val interactions = Source.fromFile(networkFile).getLines.map(line => {
      val splitted = line.split(",")
      Interaction(new Gene(splitted(0)), new Gene(splitted(1)))
    }).toSet
  new Network(interactions)
  }
}