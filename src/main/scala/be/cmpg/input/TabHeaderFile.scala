package be.cmpg.input

import java.io.InputStream

import scala.Array.canBuildFrom
import scala.io.Source

class TabHeaderFile(input: InputStream) extends Iterator[Map[String, String]]{
  val lines = Source.fromInputStream(input).getLines
  
  val names = lines.next.split("\t")
  
  def next = names.zip(lines.next.split("\t")).toMap
  def hasNext = lines.hasNext
}