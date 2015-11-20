package be.cmpg.cancer

import java.io.FileReader
import java.nio.file.Paths
import java.util.concurrent.Callable
import scala.collection.JavaConversions.asScalaBuffer
import scala.collection.JavaConversions.mutableSeqAsJavaList
import scala.collection.Set
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import au.com.bytecode.opencsv.CSVReader
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.graph.NetworkReader
import be.cmpg.walk.fungus.Fungus
import be.cmpg.utils.StatUtils
import be.cmpg.graph.Interaction
import java.io.File
import java.io.PrintWriter
import scala.collection.mutable.StringBuilder
import java.util.LinkedList
import java.util.TreeMap
import scala.collection.mutable.ListBuffer
import be.cmpg.graph.Interaction

class CancerHelper(translateGenesToEntrez: Map[String, String] = Map()) {

  def loadNetwork(networks:Seq[String]):Set[Interaction] = {
    
    var interactions = Set[Interaction]()
    val hiII = NetworkReader.fromTSV(Paths.get("networks/HI-II-14.tsv"))._1
    
    networks.foreach { n => {
      
      if (n == "HT") interactions ++= NetworkReader.fromSif(Paths.get("networks/HumanBinaryHQ_HT.txt"), 2, -1, 3)
      else if (n == "hiII14") interactions ++= hiII
      else if (n == "ppi") interactions ++= NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_ppi.txt"), 0, -1, 3)
      else if (n == "cell-map") interactions ++= NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_cell-map.sif"), 0, -1, 2)
      else if (n == "nci-nature") interactions ++= NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_nci-nature.sif"), 0, -1, 2)
      else if (n == "reactome") interactions ++= NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_reactome.sif"), 0, -1, 2)
      
      else if (n == "biogrid") interactions ++= NetworkReader.fromSif(Paths.get("networks/BIOGRID_filtered_PPI_network/PPI_network_BioGrid_HC.txt"), 0, -1, 1)
      else if (n == "kinase-substrate") interactions ++= NetworkReader.fromSif(Paths.get("networks/kinase-substrate_network/kinase-substrate_network.txt"), 0, -1, 1)
      else if (n == "encode") interactions ++= NetworkReader.fromSif(Paths.get("networks/ENCODE/regulatory.network"), 0, -1, 1)
      else {
        val networkFile = Paths.get("networks/"+n+".sif")
        if (!networkFile.toFile().exists())
          throw new RuntimeException("Network file ["+networkFile+"] was not found.")
        
        interactions ++= NetworkReader.fromSif(networkFile)
      }
    }}
    
    val hypermutatedGenes = Set("UBC", "TTN")
    interactions.
      filter( i => !hypermutatedGenes.contains(i.from.name) && !hypermutatedGenes.contains(i.to.name))
      .map(i => Interaction(i.from, i.to, if (hiII contains (i)) "hiII" else "other")).toSet
  }  
  
  def loadReferenceNetwork(referenceNetwork:String = "HT_hiII14") = {
    if (referenceNetwork == "HT_hiII14") {
      val hiII = NetworkReader.fromTSV(Paths.get("networks/HI-II-14.tsv"))._1
      val ht = NetworkReader.fromSif(Paths.get("networks/HumanBinaryHQ_HT.txt"), 2, -1, 3)
  
      val interactions = hiII ++ ht
  
      interactions.map(i => Interaction(i.from, i.to, if (hiII contains (i)) "hiII" else "ht")).toSet
      
    } else if (referenceNetwork == "MEMo") {
      val hrn1ppi_i = NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_ppi.txt"), 0, -1, 3)
      val hrn1cellmap_i = NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_cell-map.sif"), 0, -1, 2)
      val hrn1ncinature_i = NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_nci-nature.sif"), 0, -1, 2)
      val hrn1reactome_i = NetworkReader.fromSif(Paths.get("networks/MEMo/hrn1_reactome.sif"), 0, -1, 2)
  
      val interactions = hrn1ppi_i ++ hrn1cellmap_i ++ hrn1ncinature_i ++ hrn1reactome_i
  
      interactions
      
    } else if (referenceNetwork == "hiII14"){
      NetworkReader.fromTSV(Paths.get("networks/HI-II-14.tsv"))._1
      
    } else {
      throw new IllegalArgumentException("Non-implemented reference network: "+referenceNetwork)
    }
  }
  
  def loadFromDendrixFiles() = {
    val toReturn = new HashMap[(String, String), String]
    new CSVReader(new FileReader("src/test/resources/dendrix/GBM.m2"), '\t')
      //new CSVReader(new FileReader("src/test/resources/dendrix/Lung.m2"), '\t')
      //new CSVReader(new FileReader("src/test/resources/dendrix/BRCA.m2"), '\t')
      .readAll()
      .foreach(fields => {
        val sample = fields(0)
        fields.drop(1).foreach(hugoGeneName => {
          val entrezId = Gene(translateGenesToEntrez.get(hugoGeneName).getOrElse(hugoGeneName))

          toReturn.put((entrezId.name, sample), hugoGeneName)
        })
      })
    toReturn
  }

  def loadExpression(toAdd:HashMap[PolimorphismKey, Polimorphism], cnv_peaks:File, exp:File, cnv_thresholds:File) {
    println("CNV: Calculating correlation")
    val expressionR = {
      val allcnvs = readFromGeneSampleMatrix(cnv_peaks)
      val allexp = readFromGeneSampleMatrix(exp)

      val genes = allcnvs.keySet.map(_._1).toSet
      val samples = allcnvs.keySet.map(_._2).toSet

      var expressionR = new LinkedList[(String, Double)]
      for (g <- genes) {

        val X = new LinkedList[Double]
        val Y = new LinkedList[Double]

        for (s <- samples) {
          val cnv = allcnvs.get((g, s))
          val e = allexp.get((g, s))
          if (cnv.isDefined && e.isDefined) {
            X.add(cnv.get)
            Y.add(e.get)
          }
        }

        val uX = X.sum / X.size
        val uY = Y.sum / Y.size
        val sX = math.sqrt(X.map(x => math.pow(x - uX, 2)).sum)
        val sY = math.sqrt(Y.map(y => math.pow(y - uY, 2)).sum)

        val r = ((X zip Y) map { case (x, y) => (x - uX) * (y - uY) }).sum / (sX * sY)

        if (!r.isNaN() && !r.isInfinite())
          expressionR.add((g, r))
      }

      val max = expressionR.map(_._2).max

      expressionR.toList.filter(_._2 > 0.0).map { case (g, r) => (g, r / max) }.toMap // use only positive correlations abd scale so the max e
    }

    println("CNV: Reading peaks")
    
    val genesAddedToThis = new HashSet[String]
    
    var samples: LinkedList[String] = new LinkedList
    new CSVReader(new FileReader(cnv_thresholds), '\t')
      .readAll()
      .foreach(fields => {
        if (samples.isEmpty) {
          samples ++= fields.toList
        } else {
          var count = 1
          val geneName = fields(0)

          samples.foreach(sample => {
            val v = fields(count)
            val value = if (v == "NA") 0.0 else v.toDouble
            if (value.abs >= 2) {
              genesAddedToThis += geneName
              toAdd.put(PolimorphismKey(Gene(geneName), sample.replaceAll("\\.", "-")), Polimorphism(geneName, "cnv", expressionR.getOrElse(geneName, 0.0)))
            }
            count += 1
          })
        }
      })
  }
  
  def loadMaf(toAdd:HashMap[PolimorphismKey, Polimorphism], mafFile:File, dropLines:Int=2, geneNameCol:Int=0, sampleCol:Int=15, substringSample:Boolean=true, mutationTypeCol:Int=8, toExclude:List[String]=List("Silent", "Intron", "LOH")) {

    val mutations = new CSVReader(new FileReader(mafFile), '\t')
    
    for (i <- 0 until dropLines) {
      mutations.readNext()
    }
    
    var registryCount = 0
    
    var fields:Array[String] = null
    do {
      fields = mutations.readNext()
      registryCount += 1
      if (fields != null) {
        val geneName = fields(geneNameCol)
        
        var sample = fields(sampleCol)
        if (substringSample)
          sample = sample.substring(0, 15)
          
        val _type = fields(mutationTypeCol)

        if (!toExclude.contains(_type)) {
          toAdd.put(PolimorphismKey(Gene(geneName), sample.trim.replaceAll("\\.", "-")), Polimorphism(geneName, "mut", 1, _type))
        }
      }
      
      if (registryCount % 1000000 == 0) {
    	  println()
        print(registryCount+": ")        
      }
      if (registryCount % 50000 == 0)
        print('.')
      
        
    } while (fields != null)
  }
  
  def loadFromTcga2012Files() = {
	  val toReturn = new HashMap[PolimorphismKey, Polimorphism]

    println("Reading TCGA 2012: Somatic Mutations maf file")
    loadMaf(toReturn, new File("src/test/resources/tcga/BRCA2012/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf"))

    println("Reading TCGA 2012: Calculating CNV correlation with expression")

    loadExpression(toReturn, 
        new File("src/test/resources/tcga/BRCA2012/TCGA_BRCApub_cnv_peaks.txt"), 
        new File("src/test/resources/tcga/BRCA2012/BRCA.exp.547.med.txt"), 
        new File("src/test/resources/tcga/BRCA2012/TCGA_BRCApub_cnv_peaks_th.txt"))
    
    println("Reading TCGA 2012: Done.")
    toReturn.toMap

  }

  def readFromGeneSampleMatrix(file: File, dropFromHeader: Int = 0, dropFromFields: Int = 0, valueToIgnore: String = "NA") = {
    val toReturn = new HashMap[(String, String), Double]
    var samples: List[String] = List()
    val matrix = new CSVReader(new FileReader(file), '\t')

      .readAll()
      .foreach(line => {
        if (samples.isEmpty) {
          samples ++= line.drop(dropFromHeader).map(_.substring(0, 15)).toList
        } else {
          var count = 1
          val fields = line.drop(dropFromFields)
          val geneName = fields(0)
          samples.foreach(sample => {
            val v = fields(count)
            if (v != valueToIgnore && !v.isEmpty())
              toReturn.put((geneName, sample.trim.replaceAll("\\.", "-")), v.toDouble)
            count += 1
          })
        }
      })

    toReturn.view.toMap
  }

  def loadFromTcgaFiles() = {
    val toReturn = new HashMap[(String, String), Polimorphism]
    // Read somatic mutations
    new CSVReader(new FileReader("src/test/resources/tcga/BRCA/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0.somatic.maf"), '\t')
      .readAll()
      .drop(2)
      .foreach(fields => {
        val entrezId = fields(1)
        val sample = fields(15).substring(0, 15)
        toReturn.put((entrezId, sample), Polimorphism(entrezId, "snp"))

      })
    // Read hypermethylation files
    var samples = new LinkedList[String]
    new CSVReader(new FileReader("src/test/resources/tcga/BRCA/brca_tcga_hypermethylated.txt"), '\t')
      .readAll()
      .foreach(fields => {
        if (samples.isEmpty) {
          samples ++= fields.toList
        } else {
          var count = 1
          val entrezId = translateGenesToEntrez.get(fields(0)).getOrElse(fields(0))
          samples.foreach(sample => {
            val tmp = fields(count)
            if (tmp matches ("""\d+""")) {
              val value = tmp.toDouble
              if (value > 0.9) {
                toReturn.put((entrezId, sample), Polimorphism(entrezId, "hyperM"))
              }
            }
            count += 1
          })
        }
      })
    // Read hypomethylation files
    samples = new LinkedList[String]
    new CSVReader(new FileReader("src/test/resources/tcga/BRCA/brca_tcga_hypomethylated.txt"), '\t')
      .readAll()
      .foreach(fields => {
        if (samples.isEmpty) {
          samples ++= fields.toList
        } else {
          var count = 1
          val entrezId = translateGenesToEntrez.get(fields(0)).getOrElse(fields(0))
          samples.foreach(sample => {
            val tmp = fields(count)
            if (tmp matches ("""\d+""")) {
              val value = tmp.toDouble
              if (value > 0.9) {
                toReturn.put((entrezId, sample), Polimorphism(entrezId, "hypoM"))
              }
            }
            count += 1
          })
        }
      })
    // Read CNV files
    samples = new LinkedList[String]
    new CSVReader(new FileReader("src/test/resources/tcga/BRCA/Gistic_Matrix_Calculated.txt"), '\t')
      .readAll()
      .foreach(fields => {
        if (samples.isEmpty) {
          samples ++= fields.toList
        } else {
          var count = 1
          val entrezId = translateGenesToEntrez.get(fields(0)).getOrElse(fields(0))
          samples.foreach(sample => {
            val value = fields(count).toDouble
            if (value.abs >= 2) {
              toReturn.put((entrezId, sample), Polimorphism(fields(0), "exp"))
            }
            count += 1
          })
        }
      })

    toReturn
  }

  // do not use
  def loadPancancerData(workdir: String) = {

    val snpEffs = {
      val tmp = new ListBuffer[(String, Annotation)]

      // snpEffs[Chrm	Pos -> Ann(geneSymbol, impact)]
      new CSVReader(new FileReader(workdir + "/mutation.ann.canon.vcf"), '\t')
        .readAll()
        .foreach(l => {
          if (!l(0).startsWith("##") && l.length >= 7) {

            val chrm = l(0).trim
            val pos = l(1).trim
            val ann = l(7).trim.split("\\|")
            if (ann.length > 2 && (ann(1) == "downstream_gene_variant" || "MODERATE" == ann(2) || "HIGH" == ann(2))) {
            //if ("MODERATE" == ann(2) || "HIGH" == ann(2)) {
              val geneSymbol = ann(3)
              tmp += ((chrm + "_" + pos, Annotation(geneSymbol, ann(2))))
            }
          }
        })
        
      tmp.toMap
    }

    println("snpEffs: " + snpEffs.size)

    val cnvData = new File(workdir + "/cnv.tsv")
    if (!cnvData.exists()) {
      //createCnvData(workdir, cnvData)
    }
    val cnvs = new CSVReader(new FileReader(workdir + "/cnvs.txt"), '\t')
      .readAll()
      .map(l => ((l(0), l(1)), Polimorphism(l(0), "cnvs", 1.0/*l(3).toDouble*/)))

    var mutations = new CSVReader(new FileReader(workdir + "/mutations.tsv"), '\t')
    var snvs = new LinkedList[((String, String), Polimorphism)]

    var l = mutations.readNext()
    while (l != null) {
      val sample = l(0).trim
      val cancer = l(1).trim
      val chrm = l(2).trim
      val pos = l(3).trim

      val eff = snpEffs.get(chrm + "_" + pos)

      if (eff.isDefined) {
        //geneSymbol = eff.get
        snvs.add(((eff.get.geneSymbol, sample), Polimorphism(eff.get.geneSymbol, "snv", _type = cancer)))
      }
      l = mutations.readNext()
    }

    println("snvs: " + snvs.size)
    println("cnvs: " + cnvs.size)
    (snvs ++ cnvs).toMap
  }

  def createCnvData(workdir: String, cnvData: File) = {
    val printWriter = new PrintWriter(cnvData);
    println("Reading Pancancer: Calculating CNV correlation with expression")

    val expressionR = {
      val allcnvs = readFromGeneSampleMatrix(new File(workdir + "/cna.tsv"), 2, 1, "NaN")
      val allexp = readFromGeneSampleMatrix(new File(workdir + "/expression.tsv"), 2, 1, "NaN")

      val genes = allcnvs.keySet.map(_._1).toSet
      val samples = allcnvs.keySet.map(_._2).toSet

      var expressionR = new LinkedList[(String, Double)]
      for (g <- genes) {

        val X = new LinkedList[Double]
        val Y = new LinkedList[Double]

        for (s <- samples) {
          val cnv = allcnvs.get((g, s))
          val e = allexp.get((g, s))
          if (cnv.isDefined && e.isDefined) {
            X.add(cnv.get)
            Y.add(e.get)
          }
        }

        val uX = X.sum / X.size
        val uY = Y.sum / Y.size
        val sX = math.sqrt(X.map(x => math.pow(x - uX, 2)).sum)
        val sY = math.sqrt(Y.map(y => math.pow(y - uY, 2)).sum)

        val r = ((X zip Y) map { case (x, y) => (x - uX) * (y - uY) }).sum / (sX * sY)

        if (!r.isNaN() && !r.isInfinite())
          expressionR.add((g, r))
      }

      val max = expressionR.map(_._2).max

      expressionR.toList.filter(_._2 > 0.0).map { case (g, r) => (g, r / max) }.toMap // use only positive correlations abd scale so the max e
    }

    println("Reading Pancancer: CNV peaks")

    var cnas = new CSVReader(new FileReader(workdir + "/cna.tsv"), '\t')
      .readAll()

    val samples = cnas(0).drop(2).toList
    cnas
      .drop(1)
      .par
      .foreach(line => {
        var count = 1
        val fields = line.drop(1)
        val geneName = fields(0)

        samples.foreach(sample => {
          val v = fields(count)
          val value = if (v == "NaN") 0.0 else v.toDouble
          if (value.abs >= 2) {
            printWriter.println(geneName + "\t" + sample + "\t" + "cnv" + "\t" + expressionR.getOrElse(geneName, 0.0) + "\t" + 0)
            //toReturn.put((geneName, sample), Polimorphism(geneName, "cnv", expressionR.getOrElse(geneName, 0.0)))
          }
          count += 1
        })
      })

    //println("Analysed Genes: Printing")
    // gene - sample - cancerType - score - pVal
    //toReturn.foreach(e => printWriter.println (e._1._1 + "\t" + e._1._2  + "\t" + e._2.source + "\t" + e._2.score + "\t" + e._2.pVal))
    printWriter.close();
    //println("Analysed Genes: Done.")

    println("Reading Pancancer: Done!")
  }

  def loadMutationMatrixFiles(input:String):HashMap[PolimorphismKey, Polimorphism] = {
    
    
    val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]
    new CSVReader(new FileReader(input), '\t').readAll().foreach { fields =>
      val sample = fields(0)
      fields.drop(1).foreach { gene => 
        val key =   PolimorphismKey(Gene(gene), sample)
        val value = Polimorphism(gene)
        
        genePatientMatrix.put(key, value)
      }
    }    
    
    genePatientMatrix
  }
  
  def printMutationMatrixFiles(output: String, seedGenesMutations:Int, genePatientMatrix:Map[PolimorphismKey, Polimorphism]) = {
    
    val printM2 = new PrintWriter(new File(output+ ".m2"));
    genePatientMatrix.groupBy(_._1.sample).foreach { case (sample, map) =>
      if (map.size > seedGenesMutations) {
        printM2.print(sample)
        for (key <- map.keys)
          printM2.print("\t"+key.gene.name)
        
        printM2.println();
      }
    }
    genePatientMatrix.groupBy(_._1.sample)
    printM2.close();
    
    val printGlst = new PrintWriter(new File(output+ ".glst"));
    genePatientMatrix.groupBy(_._1.gene).foreach { e => 
      if (e._2.size > seedGenesMutations)
        printGlst.println(e._1.name)
    }
    printGlst.close();
  }
}

case class Annotation(val geneSymbol: String, val impact: String)
case class Polimorphism(val geneSymbol: String, val source: String = "any", val score: Double = 1, _type: String = "any")
case class PolimorphismKey(val gene:Gene, val sample:String)
case class SampleInfo(val sample:String, val subtype:String)