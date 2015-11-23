package be.cmpg.cancer

import be.cmpg.graph.NetworkReader
import be.cmpg.graph.Network
import java.nio.file.Paths
import be.cmpg.graph.Network
import scala.io.Source
import be.cmpg.graph.interaction.NetworkManager
import be.cmpg.graph.Gene
import java.util.concurrent.Callable
import be.cmpg.graph.Interaction
import be.cmpg.walk.fungus.Fungus
import be.cmpg.utils.StatUtils
import scala.collection.Set

object CancerAnalysis extends App {

  //val pheromones = List(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5);
  val pheromones = List(0.005);
  //val evaporations = List(0.99, 0.993, 0.996, 0.999, 0.9993, 0.9996, 0.9999, 1)
  val evaporations = List(0.990, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999)
  for (evaporation <- evaporations) {
    for (pheromone <- pheromones) {
      println("*********************************************")
      println("----------- Evaporation: "+evaporation)
      println("----------- Pheromones:  "+pheromone)
      println("*********************************************")
      MemoComparissonAnalysis.main(Array(2000.toString, pheromone.toString, evaporation.toString))
    }
    println();
  }

}