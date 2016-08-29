package be.cmpg.cancer

import java.io.PrintWriter
import java.io.File
import org.json.JSONObject
import org.json.JSONArray
import be.cmpg.graph.Gene
import be.cmpg.graph.interaction.NodeCostNetworkManager
import java.io.FileWriter
import java.io.FileInputStream
import java.util.HashMap
import scala.io.Source

object PrintFiles {
  // Print .m2 and .tbs mutation matrix file, the number and types of mutations per gene, the number and types of mutations per patient and a glst file containing the mutations used as input
  def printMutationMatrixFiles(output: String, seedGenesMutations: Int, genePatientMatrix: Map[PolymorphismKey, Polymorphism], weighted: Boolean, outputFolder: String) = {

    println("Printing binary (.m2) matrix")
    println("Printing samples stats");
    println("Printing as .tbs");
    val printM2 = new PrintWriter(new File(outputFolder + output + ".m2"));
    val printBySample = new PrintWriter(new File(outputFolder + output + ".bySample.stats"));
    val printTbs = new PrintWriter(new File(outputFolder + output + ".tbs"));
    genePatientMatrix.groupBy(_._1.sample).foreach {
      case (sample, map) =>
        if (map.size > seedGenesMutations) {
          printM2.print(sample)
          map.foreach {
            case (key, poli) =>
              if (weighted == true) {
                printTbs.println("column\trow\ttype\tmutation\tcnvAmp\tcnvDel\tscore")
                printM2.print("\t" + key.gene.name + "," + poli.score)
                printTbs.println(key.sample + "\t" + key.gene.name + "\t" + (if (poli.source == "mut") 1 else if (poli.source == "amp") 2 else 3) + "\t" + (if (poli.source == "mut") 1 else 0) + "\t" + (if (poli.source == "amp") 1 else 0) + "\t" + (if (poli.source == "del") 1 else 0 + "\t" + poli.score))
              } else {
                printTbs.println("column\trow\ttype\tmutation\tcnvAmp\tcnvDel")
                printM2.print("\t" + key.gene.name)
                printTbs.println(key.sample + "\t" + key.gene.name + "\t" + (if (poli.source == "mut") 1 else if (poli.source == "amp") 2 else 3) + "\t" + (if (poli.source == "mut") 1 else 0) + "\t" + (if (poli.source == "amp") 1 else 0) + "\t" + (if (poli.source == "del") 1 else 0))
              }
          }

          printM2.println();

          val bySource = map.groupBy(_._2.source).map(e => (e._1, e._2.size))
          printBySample.println(sample + "\t" + map.size + "\t" + bySource.mkString(","))
        }
    }
    printBySample.close();
    printM2.close();
    printTbs.close();

    println("Printing genes and genes stats");

    {
      val printByGene = new PrintWriter(new File(outputFolder + output + ".byGene.stats"));
      val printGlst = new PrintWriter(new File(outputFolder + output + ".glst"));
      genePatientMatrix.groupBy(_._1.gene).foreach { e =>
        if (e._2.size >= seedGenesMutations) {
          val bySource = e._2.groupBy(_._2.source).map(e => (e._1, e._2.size))
          printGlst.println(e._1.name)
          printByGene.println(e._1.name + "\t" + e._2.size + "\t" + bySource.mkString(","))
        }
      }
      printGlst.close();
      printByGene.close();
    }

  }
  
  // This function writes the html network file as well as the _edges.tsv and the _nodes.tsv and the _pattern.tsv files
def printPattern(outputPrefix: String, genes: List[Gene], networkManager: NodeCostNetworkManager, genePatientMatrix: Map[PolymorphismKey, Polymorphism], outputFolder:String, inputFolder:String, additional: Map[Gene, Map[String, Any]] = Map(), otherPositiveGeneSetLists:Set[Gene]=Set()) = {
                                        
    val results = new JSONObject

    val edges = new JSONArray

    val additionalNodeInfo = if (additional.isEmpty && !genePatientMatrix.isEmpty) {
      val mutationsPerGene = genePatientMatrix.groupBy(_._1.gene)

      genes.map { g =>
        val info = mutationsPerGene.get(g)
        if (info.isDefined)
          g -> Map[String, Any]("pvalue" -> math.sqrt(info.get.size), "origin" -> "coding")
        else
          g -> Map[String, Any]("pvalue" -> 0.0, "origin" -> "unkown")
      }.toMap[Gene, Map[String, Any]]
    } else {
       additional 
    }

    val out_edges = new FileWriter(outputFolder+outputPrefix + "_edges.tsv")
    out_edges.write("GeneA\tGeneB\tMEScore\n")
    networkManager.getAllInteractions
      .filter { i => i.to != i.from && genes.contains(i.from) && genes.contains(i.to) }
      .foreach { i =>
        {
          val score = { val score = networkManager.scoreSubnetwork(Set(i)); if (score.isNaN) 0.0 else score }
          
          val edgeInfo = new JSONObject()
            .accumulate("source", i.from.name)
            .accumulate("target", i.to.name)
            .accumulate("score", score)
            .accumulate("evidence",i.evidence.mkString(","))
          edges.put(edgeInfo)
          out_edges.write(i.from.name + "\t" + i.to.name + "\t" + score + "\n")
        }
      }

    out_edges.close()
    
    def getKnownCancerLists(gene:Gene) = {
      var toReturn = 0;
      
      if (otherPositiveGeneSetLists.contains(gene)) toReturn += 1
      if (DataLoader.cgcGenes(inputFolder).contains(gene)) toReturn += 1
      if (DataLoader.ncgGenes(inputFolder).contains(gene)) toReturn += 1
      
      toReturn
    }  
    
    val truePositiveGeneList = DataLoader.cgcGenes(inputFolder) ++ DataLoader.ncgGenes(inputFolder) ++ otherPositiveGeneSetLists;
    
    val nodes = new JSONArray

    val out_nodes = new FileWriter(outputFolder+outputPrefix + "_nodes.tsv")
    out_nodes.write("GeneSymbol\tPosteriorP\tSelected\tConvergenceIter\tKnownCancerGene\t\tBestSSN\n")
    println("GeneSymbol\tOtherPositive\tCGC\tNCG")
    genes
      .foreach(g => {
        
        println(g.name+"\t"+otherPositiveGeneSetLists.contains(g)+"\t"+DataLoader.cgcGenes(inputFolder).contains(g)+"\t"+DataLoader.ncgGenes(inputFolder).contains(g))
        
        val geneInfo = new JSONObject()
          .accumulate("name", g.name)
          .accumulate("selected", genes.contains(g))
          .accumulate("knownCancerGene",("["+(DataLoader.cgcGenes(inputFolder).contains(g),otherPositiveGeneSetLists.contains(g),DataLoader.ncgGenes(inputFolder).contains(g)).toString+"]").replace("(","").replace(")","").replace("\"",""))
          .accumulate("inTPLists", getKnownCancerLists(g))
          
        val additionalInfo = additionalNodeInfo.get(g)
        if (additionalInfo.isDefined) {
          additionalInfo.get.foreach { e =>
            geneInfo.accumulate(e._1, e._2)
          }
        }
        
        nodes.put(geneInfo)
        try {
        	val node = networkManager.getNetwork().getNode(g)
        	out_nodes.write(g.name + "\t" + node.posteriorProbability+ "\t" + genes.contains(g) + "\t" + node.convergenceIteration + "\t" + truePositiveGeneList.contains(g) + "\t" + node.bestSubnetwork._1+"\t"+node.bestSubnetwork._2+"\n")
        } catch {
          case t: Throwable => out_nodes.write(g.name + "\t" + "NA" + "\t" + genes.contains(g) + "\t" + "NA" + "\t" + truePositiveGeneList.contains(g) + "\t" + "NA" + "\n")
        }
      })

    out_nodes.close()

    results.put("edges", edges)
    results.put("nodes", nodes)

    val html = Source.fromInputStream(new FileInputStream(inputFolder+"NetworkOutput.html")).mkString

    val out_json = new FileWriter(outputFolder+outputPrefix + "_network.html")
    out_json.write(html.replace("$${graph}$$", results.toString.replace("\"[","[").replace("]\"","]")))

    out_json.close()

    val all_samples = genePatientMatrix.map(_._1.sample).toSet

    val geneMutScoreCounts = genes.map(gene => (gene, {
      var mutationCounts = 0
      for (sample <- all_samples) {
        // Location not yet implemented here !
        if (genePatientMatrix contains (PolymorphismKey(gene, sample,0))) {
          mutationCounts += 1
        }
      }
      var mutualExScore = 0.0
      for (sample <- all_samples) {
        // Location not yet implemented here !
        if (genePatientMatrix contains (PolymorphismKey(gene, sample,0))) {

          var nonMutualGenes = 0.0
          for (other <- genes) {
            // Location not yet implemented here !
            if (genePatientMatrix.contains(PolymorphismKey(other, sample,0))) {
              nonMutualGenes += 1
            }
          }

          mutualExScore += 1.0 / nonMutualGenes
        }
      }

      val rank = genes.indexOf(gene).toDouble
      (mutualExScore, mutationCounts, rank)
    })).toList.sortBy(_._2._3)

    val toPrint = new HashMap[String, List[String]]

    val out = new FileWriter(outputFolder+outputPrefix + "_pattern.tsv")

    val header = {
      out.append("Gene")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._1.name))
      out.append("\n").append("Rank")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._2._3.toString()))
      out.append("\n").append("GeneName")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._1.name))
      out.append("\n").append("MutualExclScore")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._2._1.toString()))
      out.append("\n").append("MutQty")
      geneMutScoreCounts.foreach(g => out.append("\t").append(g._2._2.toString()))
    }

    val samplesStr = all_samples.map(sample => (sample, {
      val sb = new StringBuilder
      // Location not yet implemented here !
      geneMutScoreCounts.foreach(gMsc => sb.append("\t").append(genePatientMatrix.getOrElse(PolymorphismKey(gMsc._1, sample,0), Polymorphism("Nothing", score = 9)).score.toShort))
      sb.toString
    })).toList.sortBy(_._2)

    val mutualExclusivityPattern = {
      samplesStr.foreach(s => out.append("\n").append(s._1).append(s._2.replaceAll("9", "")))
    }
  }
  
}
