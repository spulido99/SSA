package be.cmpg.graph

import java.nio.file.Path
import scala.io.Source
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.JavaConversions._
import au.com.bytecode.opencsv.CSVParser

object NetworkReader {

  def fromSif(file: Path, fromIdx:Int = 0, typeIdx:Int = 1, toIdx:Int = 2): Set[Interaction] = {
    new CSVReader(new FileReader(file.toFile()), '\t')
    .readAll()
    .map(fields => {
      val from = Gene(fields(fromIdx))
      val iType = if (typeIdx < 0) "undirected" else fields(typeIdx)
      val to = Gene(fields(toIdx))
      Interaction(from, to, iType)
    }).toSet
  }
  
  def fromFile(file: Path): Set[Interaction] = {

    println("*********************** Reading network from File")

    //CSVReader(Reader reader, char separator, char quotechar, char escape, int line, boolean strictQuotes, boolean ignoreLeadingWhiteSpace)
    val reader = new CSVReader(new FileReader(file.toFile()), '\t');

    val interactions = new HashSet[Interaction]()

    var first = true
    for (fields <- reader.readAll()) {
      if (!fields(0).startsWith("%")){
        val direction = if (fields(2) == "pp" || fields(2) == "met_undirected") "undirected" else "directed"
        val regulatory = if (Set("pd", "srna", "sigma").contains(fields(2))) "regulatory" else "non-regulatory"
        val prob = if (fields.size > 3) fields(3).toDouble else 1.0
        Interaction(Gene(fields(0)), Gene(fields(1)), fields(2), direction, regulatory, prob)
      }
    }

    interactions.toSet
    
    Source.fromURI(file.toUri())
      .getLines()
      .filter(!_.startsWith("%"))
      .map(line => {
        val splitted = line.split("\t")
        val direction = if (splitted(2) == "pp" || splitted(2) == "met_undirected") "undirected" else "directed"
        val regulatory = if (Set("pd", "srna", "sigma").contains(splitted(2))) "regulatory" else "non-regulatory"
        val prob = if (splitted.size > 3) splitted(3).toDouble else 1.0
        Interaction(Gene(splitted(0)), Gene(splitted(1)), splitted(2), direction, regulatory, prob)
      })
      .toSet
  }

  def fromTSV(file: Path, typ:Option[String] = None) : (Set[Interaction], Map[String, String]) = {

    val interactions = new HashSet[Interaction]()
    interactions.useSizeMap(true)
    
    val translateGenesToEntrez = new HashMap[String, String]
    
    val reader = new CSVReader(new FileReader(file.toFile()), '\t');
    
    var first = true
    var edgeCount = 0
    for (fields <- reader.readAll()) {
      if (!first) {
        edgeCount += 1
        
        interactions += Interaction(Gene(fields(1).trim), Gene(fields(3).trim), if (typ.isDefined) typ.get else "any")
        
        translateGenesToEntrez.put(fields(0).trim, fields(1).trim)
        translateGenesToEntrez.put(fields(2).trim, fields(3).trim)
        
      } else {
        first = false
      }
    }
    (interactions.toSet, translateGenesToEntrez.toMap)
  }
  
  def fromCytoScapeFiles(edges: Path, nodes: Path): (Set[Interaction], Map[String, String]) = {

    val interactions = new HashSet[Interaction]()
    interactions.sizeHint(160000)
    interactions.useSizeMap(true)

    val translateGenesToEntrez = new HashMap[String, String]

    println("*********************** Reading from cytoscape files [Edges]")
    val edgesReader = new CSVReader(new FileReader(edges.toFile()));
    //val myEntries = reader.readAll()
    //val edgesReader = CSVReader.open()

    var first = true
    var edgeCount = 0
    for (fields <- edgesReader.readAll()) {
      if (!first) {
        edgeCount += 1
        val interactionStr = fields(2).split(" ")
        interactions += Interaction(Gene(interactionStr(0).trim), Gene(interactionStr(2).trim))
        if (edgeCount % 500 == 0)
          print(".")
        if (edgeCount % 20000 == 0)
          println()
      } else {
        first = false
      }
    }
    println()
    println(interactions.size + " edges")

    println("*********************** Reading from cytoscape files [Nodes]")
    //val nodesReader = CSVReader.open(nodes.toFile())
    val nodesReader = new CSVReader(new FileReader(nodes.toFile()))

    first = true

    var nodeCount = 0

    for (fields <- nodesReader.readAll()) {
      if (!first) {
        nodeCount += 1
        val entrezId = fields(1)
        val aliases = fields(4).split("\n")
        for (alias <- aliases) {
            
          if (translateGenesToEntrez.contains(alias)) {
            //println("Gene alias repeated: "+alias+". Alreade asociated with entrezId "+translateGenesToEntrez(alias))
          } else {
            translateGenesToEntrez.put(alias, entrezId)
          }
        }
      } else {
        first = false
      }
    }
    println(nodeCount + " nodes")

    (interactions.toSet, translateGenesToEntrez.toMap)
  }
}