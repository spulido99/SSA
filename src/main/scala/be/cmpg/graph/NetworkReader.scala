package be.cmpg.graph

import java.nio.file.Path
import scala.io.Source
import scala.collection.mutable.HashSet
import scala.collection.mutable.HashMap
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.JavaConversions._
import au.com.bytecode.opencsv.CSVParser
import java.io.InputStream
import java.net.URL
import java.io.BufferedReader
import java.io.Reader
import java.io.InputStreamReader
import java.nio.file.Paths
import java.io.File

object NetworkReader {

  def getReader(fileName:String) = {
    val file = new File(fileName)
    if (file.exists()) 
       new FileReader(file) 
    else 
    	 new InputStreamReader(getClass.getResourceAsStream("/"+fileName))
  }
  
  def fromSif(file: String, fromIdx:Int, typeIdx:Int, toIdx:Int, evidenceId:Int= -1): Set[Interaction] = {
    
    new CSVReader(getReader(file), '\t')
    .readAll()
    .drop(1)
    .map(fields => {
      val from = Gene(fields(fromIdx))
      val iType = if (typeIdx < 0) "undirected" else fields(typeIdx)
      val to = Gene(fields(toIdx))
      val evidence = if(evidenceId== -1){Set("Not available")} else{fields(evidenceId).split(",").toSet}
//      Interaction(from, to, iType)
      Interaction(from, to, iType,evidence=evidence)
    }).toSet
  }
  
  def fromFile(file:String): Set[Interaction] = {

    println("*********************** Reading network from File")

    //CSVReader(Reader reader, char separator, char quotechar, char escape, int line, boolean strictQuotes, boolean ignoreLeadingWhiteSpace)
    val reader = new CSVReader(getReader(file), '\t');

    val interactions = new HashSet[Interaction]()

    var first = true
    for (fields <- reader.readAll()) {
      if (!fields(0).startsWith("%")){
        val direction = if (fields(2) == "pp" || fields(2) == "met_undirected") "undirected" else "directed"
        val regulatory = if (Set("pd", "srna", "sigma").contains(fields(2))) "regulatory" else "non-regulatory"
        val prob = if (fields.size > 3) fields(3).toDouble else 1.0
//        Interaction(Gene(fields(0)), Gene(fields(1)), fields(2), direction, regulatory, prob)
        Interaction(Gene(fields(0)), Gene(fields(1)), fields(2), direction, regulatory)
      }
    }

    interactions.toSet
    
    Source.fromURI(Paths.get(file).toUri())
      .getLines()
      .filter(!_.startsWith("%"))
      .map(line => {
        val splitted = line.split("\t")
        val direction = if (splitted(2) == "pp" || splitted(2) == "met_undirected") "undirected" else "directed"
        val regulatory = if (Set("pd", "srna", "sigma").contains(splitted(2))) "regulatory" else "non-regulatory"
        val prob = if (splitted.size > 3) splitted(3).toDouble else 1.0
 //       Interaction(Gene(splitted(0)), Gene(splitted(1)), splitted(2), direction, regulatory, prob)
        Interaction(Gene(splitted(0)), Gene(splitted(1)), splitted(2), direction, regulatory)
      })
      .toSet
  }

  def fromTSV(file:String, typ:Option[String] = None) : (Set[Interaction], Map[String, String]) = {

    val interactions = new HashSet[Interaction]()
    interactions.useSizeMap(true)
    
    val translateGenesToEntrez = new HashMap[String, String]    
    
    val reader = new CSVReader(getReader(file), '\t');
    
    var first = true
    var edgeCount = 0
    for (fields <- reader.readAll()) {
      if (!first) {
        edgeCount += 1
    
    val interactionEvidence = if (!(fields(4)=="")){fields(4).split(",").toSet} else {Set("Not available")}
       
        interactions += Interaction(Gene(fields(1).trim), Gene(fields(3).trim), if (typ.isDefined) typ.get else "any",evidence=interactionEvidence)
        
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