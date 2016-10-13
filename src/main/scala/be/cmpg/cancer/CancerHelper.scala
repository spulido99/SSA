package be.cmpg.cancer

import java.io.FileReader
import java.nio.file.Paths
import java.util.concurrent.Callable
import scala.collection.JavaConversions.asScalaBuffer
import scala.collection.JavaConversions.mutableSeqAsJavaList
import scala.collection.Set
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import scala.collection.mutable.ArrayBuffer
import au.com.bytecode.opencsv.CSVReader
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.graph.Network
import be.cmpg.graph.NetworkReader
import be.cmpg.walk.fungus.Fungus
import be.cmpg.utils.StatUtils
import java.io.File
import java.io.PrintWriter
import scala.collection.mutable.StringBuilder
import java.util.LinkedList
import java.util.TreeMap
import scala.collection.mutable.ListBuffer
import scopt.OptionParser
import scala.util.Random
import scala.util.Try
import scala.io.Source
import java.io.InputStreamReader

class CancerHelper(translateGenesToEntrez: Map[String, String] = Map()) {

  val random = new Random(System.nanoTime())
  
  lazy val cgc = new CSVReader(new InputStreamReader(getClass.getResourceAsStream("/cancergenelists/Census_allFri20150102.csv"))) // CGC: Cancer Gene Census
                                .readAll()
                                .drop(1)
                                .map(_(0))
                                .map(Gene(_))
                                .toSet    
    
  lazy val ncg = new CSVReader(new InputStreamReader(getClass.getResourceAsStream("/cancergenelists/ncg.txt")), '\t') // NCG 4.0: Network of Cancer Genes
                                .readAll()
                                .drop(1)
                                //.filter { x => x(25) == "FALSE" }
                                .map(_(1))
                                .map(Gene(_))
                                .toSet
                                
  lazy val zackGenes = Source.fromInputStream(getClass.getResourceAsStream("/cancergenelists/ZackEtAl_CNAGenes.txt")).getLines.toSet

  def loadNetwork(networks: Seq[String]): Set[Interaction] = {

    var interactions = Set[Interaction]()
    val hiII = NetworkReader.fromTSV("networks/HI-II-14.tsv")._1

    networks.foreach { n =>
      {

        if (n == "HT") interactions ++= NetworkReader.fromSif("networks/HumanBinaryHQ_HT.txt", 2, -1, 3)
        else if (n == "hiII14") interactions ++= hiII
        else if (n == "ppi") interactions ++= NetworkReader.fromSif("networks/MEMo/hrn1_ppi.txt", 0, -1, 3)
        else if (n == "cell-map") interactions ++= NetworkReader.fromSif("networks/MEMo/hrn1_cell-map.sif", 0, -1, 2)
        else if (n == "nci-nature") interactions ++= NetworkReader.fromSif("networks/MEMo/hrn1_nci-nature.sif", 0, -1, 2)
        else if (n == "reactome") interactions ++= NetworkReader.fromSif("networks/MEMo/hrn1_reactome.sif", 0, -1, 2)

        else if (n == "biogrid") interactions ++= NetworkReader.fromSif("networks/BIOGRID_filtered_PPI_network/PPI_network_BioGrid_HC.txt", 0, -1, 1)
        else if (n == "kinase-substrate") interactions ++= NetworkReader.fromSif("networks/kinase-substrate_network/kinase-substrate_network.txt", 0, -1, 1)
        else if (n == "encode") interactions ++= NetworkReader.fromSif("networks/ENCODE/regulatory.network", 0, -1, 1)
        else {
          val networkFile = "networks/" + n + ".sif"
          if (!new File(networkFile).exists())
            throw new RuntimeException("Network file [" + networkFile + "] was not found.")

          interactions ++= NetworkReader.fromSif(networkFile)
        }
      }
    }

    val hypermutatedGenes: Set[String] = Set();
    interactions.
      filter(i => !hypermutatedGenes.contains(i.from.name) && !hypermutatedGenes.contains(i.to.name))
      .map(i => Interaction(i.from, i.to, if (hiII contains (i)) "hiII" else "other")).toSet
  }

  def loadReferenceNetwork(referenceNetwork: String = "HT_hiII14") = {
    if (referenceNetwork == "HT_hiII14") {
      val hiII = NetworkReader.fromTSV("networks/HI-II-14.tsv")._1
      val ht = NetworkReader.fromSif("networks/HumanBinaryHQ_HT.txt", 2, -1, 3)

      val interactions = hiII ++ ht

      interactions.map(i => Interaction(i.from, i.to, if (hiII contains (i)) "hiII" else "ht")).toSet

    } else if (referenceNetwork == "MEMo") {

      val hrn1ppi_i = NetworkReader.fromSif("networks/MEMo/hrn1_ppi.txt", 0, -1, 3)
      val hrn1cellmap_i = NetworkReader.fromSif("networks/MEMo/hrn1_cell-map.sif", 0, -1, 2)
      val hrn1ncinature_i = NetworkReader.fromSif("networks/MEMo/hrn1_nci-nature.sif", 0, -1, 2)
      val hrn1reactome_i = NetworkReader.fromSif("networks/MEMo/hrn1_reactome.sif", 0, -1, 2)

      val interactions = hrn1ppi_i ++ hrn1cellmap_i ++ hrn1ncinature_i ++ hrn1reactome_i

      interactions

    } else if (referenceNetwork == "hiII14") {
      NetworkReader.fromTSV("networks/HI-II-14.tsv")._1

    } else {
      throw new IllegalArgumentException("Non-implemented reference network: " + referenceNetwork)
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

  def loadExpression(cnvdata: Map[String, File], acceptedCorrelationQVal: Double = 0.05, gisticThreshold:Int=1, maxQtyCopyNumber:Int = 500):Map[PolimorphismKey, Polimorphism] = {

    val gistic = cnvdata.get("gistic")
    val peaks = cnvdata.get("cnv_peaks")
    val expression = cnvdata.get("exp")
    val correlation = cnvdata.get("corr")

    if (gistic.isEmpty)
      throw new RuntimeException("If using CNV gistic folder should exist.")

    val expressionR: Map[String, Double] = if (correlation.isDefined) {
      println("CNV: loading correlation")

      new CSVReader(new FileReader(correlation.get), '\t')
        .readAll()
        .filter { fields => Try { fields(3).toDouble }.getOrElse(1.0) <= acceptedCorrelationQVal } // Column (field) 3 is the qvalue
        .map { fields =>
          val gene = fields(0)
          val corr = Try { fields(1).toDouble }.getOrElse(0.0)
          (gene, corr)
        }.toMap

    } else if (expression.isDefined && peaks.isDefined) {

      println("CNV: Calculating correlation")

      val allcnvs = readFromGeneSampleMatrix(peaks.get)
      val allexp = readFromGeneSampleMatrix(expression.get)

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

    } else {
      throw new RuntimeException("No Correlation or files to calculate correlation available.")
    }

    println("CNV: Reading peaks")


    //var samples: List[String] = List()
    val nonsamples = List("Gene Symbol", "Locus ID", "Cytoband")
    
    val gisticAmpGenes = loadGisticConfidenceGenes(gistic.get + "/amp_genes.conf_99.txt")
    
    val gisticDelGenes = loadGisticConfidenceGenes(gistic.get + "/del_genes.conf_99.txt")
    
    val confGenes = (gisticAmpGenes ++ gisticDelGenes).toSet intersect (zackGenes ++ cgc ++ ncg)
    
    val gisticOutput = Source.fromFile(gistic.get + "/all_thresholded.by_genes.txt", "latin1").getLines.map { _.split("\t") }
    
    val samples = gisticOutput
            .next()
            .toList
            .filterNot { nonsamples contains _ } // remove those that are not samples
            .map { s => // change the names to sample size if they are TCGA samples
              var sample = s.replaceAll("\\.", "-")
              if (sample.startsWith("TCGA-") && sample.size > 15) {
                sample = sample.substring(0, 15)
              }
              sample
            }
    
    val all = gisticOutput.map { fields => 
        var count = fields.size - samples.size - 1
        val geneName = fields(0)
        samples.map {sample => 
      	  count += 1
          if (!confGenes.contains(geneName)) {
            None
          } else if (nonsamples.contains(sample)) {
          	None
          } else {
          	val value = Try { fields(count).toInt }.getOrElse(0)
      			val geneExpressionCorr = expressionR.get(geneName)
      			if (value.abs < gisticThreshold || geneExpressionCorr.isEmpty) {
      				None
            } else {
          	  Some(PolimorphismKey(Gene(geneName), sample), Polimorphism(geneName, if (value > 0) "amp" else "del", geneExpressionCorr.get))
            }
          }
        }
    }.flatten.flatten.toList
    
    val hypermutators = all.groupBy(_._1.sample).filter( _._2.size > maxQtyCopyNumber ).keys.toSet
    
    all.filter { e => !hypermutators.contains(e._1.sample) }.toMap

    
    /*
    new CSVReader(new FileReader(thresholds.get), '\t')
      .readAll()
      .foreach(fields => {
        if (samples.isEmpty) {
          samples = fields
            .toList
            .filterNot { nonsamples contains _ } // remove those that are not samples
            .map { s => // change the names to sample size if they are TCGA samples
              var sample = s.replaceAll("\\.", "-")
              if (sample.startsWith("TCGA-") && sample.size > 15) {
                sample = sample.substring(0, 15)
              }
              sample
            }
        } else {
          var count = fields.size - samples.size
          val geneName = fields(0)

          samples.foreach(sample => {
            if (!nonsamples.contains(sample)) {
              val value = Try { fields(count).toDouble }.getOrElse(0.0)
              val geneExpressionCorr = expressionR.get(geneName)
              if (value.abs >= 2 && geneExpressionCorr.isDefined) {
                genesAddedToThis += geneName
                toAdd.put(PolimorphismKey(Gene(geneName), sample), Polimorphism(geneName, "cnv", geneExpressionCorr.get))
              }
            }
            count += 1
          })
        }
      })*/
  }
  
  def loadGisticConfidenceGenes(fileName:String) = {
    
    Source.fromFile(fileName, "latin1")
                        .getLines
                        .drop(4)
                        .map { x => x.split("\t") }
                        .flatten
                        .filter { x => x.length() > 0 }
  }

  def loadMaf(mafFile: File, dropLines: Int = 2, geneNameCol: Int = 0, sampleCol: Int = 15, substringSample: Boolean = true, mutationTypeCol: Int = 8, toExclude: List[String] = List("Silent", "Intron", "LOH"), maxQtyMutations:Int=500):Map[PolimorphismKey, Polimorphism] = {

  val all = Source.fromFile(mafFile, "latin1")
        .getLines
        .drop(dropLines)
        .map { _.split("\t") }
        .map { fields =>

          val geneName = fields(geneNameCol)
          var sample = fields(sampleCol)
          if (substringSample)
            sample = sample.substring(0, 15)
  
          val _type = fields(mutationTypeCol)
  
          if (toExclude.contains(_type)) {
        	  None
          } else {
        	  Some(PolimorphismKey(Gene(geneName), sample.trim.replaceAll("\\.", "-")), Polimorphism(geneName, "mut", 1, _type))
          }
        }
        .flatten
        .toList

  val hypermutators = all.groupBy(_._1.sample).filter( _._2.size > maxQtyMutations ).keys.toSet
  
  all.filter { e => !hypermutators.contains(e._1.sample) }.toMap
  
  /*
    
    val mutations = new CSVReader(new FileReader(mafFile), '\t')

    
    for (i <- 0 until dropLines) {
      mutations.readNext()
    }

    var registryCount = 0

    var fields: Array[String] = null
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
        print(registryCount + ": ")
      }
      if (registryCount % 50000 == 0)
        print('.')

    } while (fields != null)*/
  }

  def loadFromTcga2012Files() = {

    println("Reading TCGA 2012: Somatic Mutations maf file")
    val toReturnM = loadMaf(new File("src/test/resources/tcga/BRCA2012/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf"))

    println("Reading TCGA 2012: Calculating CNV correlation with expression")

    val toReturnE = loadExpression(Map(
      "cnv_peaks" -> new File("src/test/resources/tcga/BRCA2012/TCGA_BRCApub_cnv_peaks.txt"),
      "exp" -> new File("src/test/resources/tcga/BRCA2012/BRCA.exp.547.med.txt"),
      "cnv_peaks" -> new File("src/test/resources/tcga/BRCA2012/TCGA_BRCApub_cnv_peaks_th.txt")))

    println("Reading TCGA 2012: Done.")
    toReturnM ++ toReturnE
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
      .map(l => ((l(0), l(1)), Polimorphism(l(0), "cnvs", 1.0 /*l(3).toDouble*/ )))

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

  def loadMutationMatrixFiles(input: String, hyperMutators:Int=400): HashMap[PolimorphismKey, Polimorphism] = {

    val genePatientMatrix = new HashMap[PolimorphismKey, Polimorphism]
    new CSVReader(new FileReader(input), '\t').readAll().foreach { fields =>
      
      if (fields.length < hyperMutators) {
        val sample = fields(0)
        fields.drop(1).foreach { gene =>
          val key = PolimorphismKey(Gene(gene), sample)
          val value = Polimorphism(gene)
    
          genePatientMatrix.put(key, value)
        }
      }
    }

    genePatientMatrix
  }

  def printMutationMatrixFiles(output: String, seedGenesMutations: Int, genePatientMatrix: Map[PolimorphismKey, Polimorphism]) = {

    println("Printing binary (.m2) matrix")
    println("Printing samples stats");
    println("Printing as .tbs");
    val printM2 = new PrintWriter(new File(output + ".m2"));
    val printBySample = new PrintWriter(new File(output + ".bySample.stats"));
    val printTbs = new PrintWriter(new File(output + ".tbs"));
    printTbs.println("column\trow\ttype\tmutation\tcnvAmp\tcnvDel")
    genePatientMatrix.groupBy(_._1.sample).foreach {
      case (sample, map) =>
        if (map.size > seedGenesMutations) {
          printM2.print(sample)
          map.foreach{ case (key, poli) => 
            printM2.print("\t" + key.gene.name)
            printTbs.println(key.sample + "\t" + key.gene.name + "\t" + (if (poli.source=="mut") 1 else if (poli.source=="amp") 2 else 3) + "\t" + (if (poli.source=="mut") 1 else 0) + "\t" +(if (poli.source=="amp") 1 else 0)+ "\t" +(if (poli.source=="del") 1 else 0))
          }

          printM2.println();
        
          val bySource = map.groupBy(_._2.source).map(e => (e._1, e._2.size))
          printBySample.println(sample + "\t" + map.size + "\t"+bySource.mkString(","))
        }
    }
    printBySample.close();
    printM2.close();
    printTbs.close();

    println("Printing genes and genes stats");
    
    {
      val printByGene = new PrintWriter(new File(output + ".byGene.stats"));
      val printGlst = new PrintWriter(new File(output + ".glst"));
      genePatientMatrix.groupBy(_._1.gene).foreach { e =>
        if (e._2.size > seedGenesMutations) {
          val bySource = e._2.groupBy(_._2.source).map(e => (e._1, e._2.size))
          printGlst.println(e._1.name)
          printByGene.println(e._1.name + "\t" + e._2.size+ "\t" + bySource.mkString(","))
        }
      }
      printGlst.close();
      printByGene.close();
    }
    
  }

  def getBasicArgParser(name: String): OptionParser[Config] = {

    new scopt.OptionParser[Config](name) {

      head("Mutual exclusivity analysis by Small Subnetwork Analysis", "0.1")
      opt[Int]('i', "iterations") action { (x, c) =>
        c.copy(iterations = x)
      } text ("iterations is the number of iterations (default: 50000)")

      opt[Double]('r', "reinforcement") action { (x, c) =>
        c.copy(reinforcement = x)
      } text ("The reinforncement for succesful learning (default: 0.0005)")

      opt[Double]('f', "forgetfulness") action { (x, c) =>
        c.copy(forgetfulness = x)
      } text ("The forgetfulness after each iteration (default: 0.9996)")

      opt[Int]('s', "seedGenesMutations") action { (x, c) =>
        c.copy(seedGenesMutations = x)
      } text ("The number of mutated samples required in a gene to be a seed (default: 1)")

      opt[String]('o', "outputPrefix") action { (x, c) =>
        c.copy(outputPrefix = x)
      } text ("The prefix for the analisis output files (default: ME)")

      opt[Seq[String]]('n', "refNetwork") action { (x, c) =>
        c.copy(refNetwork = x)
      } text ("The reference network (HT, hiII14, ppi, cell-map, nci-nature, reactome or/and files. HT and hiII14 used by default.)")

      /*    
    opt[Int]('g', "maxOutputGenes") action { (x, c) =>
      c.copy(maxOutputGenes = x)
    } text ("Jump by ")
    */

      opt[Boolean]('u', "useRank") action { (x, c) =>
        c.copy(useRank = x)
      } text ("Rank the gene scores instead of the scaled values (default: true)")

      opt[Boolean]('h', "statistical") action { (x, c) =>
        c.copy(statistical = x)
      } text ("Use an hypergeometric test to search for mutual exclusivity. Will not find low mutations genes. (default: false)")

      opt[Int]('g', "minMutPerGene") action { (x, c) =>
        c.copy(minMutPerGene = x)
      } text ("Ignore genes mutated in less than 'g' samples (default: 1)")

      opt[String]('m', "input_matrix") required () action { (x, c) =>
        var inputMatrix = x
        if (!x.endsWith(".m2"))
          inputMatrix = x + ".m2"
        c.copy(input = Some(inputMatrix))
      } text ("Mutation matrix .m2 file (created with option ME_input )")

      opt[Int]('p', "outputGenes") action { (x, c) =>
        c.copy(outputGenes = x)
      } text ("Number of nodes to be published to the html network output (default: 100).")

      opt[Int]("processors") action { (x, c) =>
        c.copy(processors = x)
      } text ("Number parallel processors to use (default: half of available processors).")

      opt[File]("subtypeFile") action { (x, c) =>
        c.copy(subtypeFile = Some(x))
      } text ("Samples per subtype (Tab Delimited File with 2 columns: Sample Subtype")

      opt[String]("subtype") action { (x, c) =>
        c.copy(subtype = Some(x))
      } text ("subtype to analyse")
      
      opt[Seq[String]]("excludeGenes") action { (x, c) =>
        c.copy(excludeGenes = Some(x.map { Gene(_) }))
      } text ("Genes to remove from the network analysis (e.g. --excludeGenes TP53,MYC)")

      opt[Seq[String]]("debug") action { (x, c) =>
        c.copy(debug = Some(x))
      } text ("Debug mode. List of genes that want to be observed (e.g. -d TP53,MYC)")
      
      opt[Boolean]("randomizeData") action { (x, c) =>
        c.copy(randomizeData = x)
      } text ("It will randomize the input data by genenames. Useful to see differences in small subnetworks score distributions. (default: false)")
      
      

      /*
      opt[Boolean]("calculatePValue") action { (x, c) =>
      c.copy(calculatePValue = x)
      } text ("Calculates p-values for the selected genes. Create a new output file with an experimental p-value per gene. (argument \"--outputGenes\" is required)")
      
      opt[Int]("pvalIterations") action { (x, c) =>
        c.copy(pvalIterations = x)
      } text ("Different number of iterations can be used to speed up p-value calculations. (Only used for p-value calculations)")

      opt[Double]("pvalReinforcement") action { (x, c) =>
        c.copy(pvalReinforcement = x)
      } text ("Different values of reinforcement can be used to speed up p-value calculations. (Only used for p-value calculations)")

      opt[Double]("pvalForgetfulness") action { (x, c) =>
        c.copy(pvalForgetfulness = x)
      } text ("Different values of forgetfulness can be used to speed up p-value calculations. (Only used for p-value calculations)")
      */
      help("help")
    }
  }

  def loadData(config: Config) = {
    var interactions = this.loadNetwork(config.refNetwork)
    if (config.excludeGenes.isDefined) {
    	interactions = interactions.filter { i => !config.excludeGenes.get.contains(i.from) && !config.excludeGenes.get.contains(i.to) }
    }
    
    val network = new Network(interactions)

    println("Interactions: " + interactions.size)

    val genePatientMatrix = {

      println("Loading mutation matrix file...")
      val genePatientMatrix = this.loadMutationMatrixFiles(config.input.get)

      if (config.subtype.isDefined && config.subtypeFile.isDefined) {

        //subtypeFile = TCGA_BRCA_clin_PAM50.txt
        val samplesInfo = new CSVReader(new FileReader(config.subtypeFile.get), '\t')
          .readAll()
          .map { fields => SampleInfo(fields(0).trim.replaceAll("\\.", "-"), fields(fields.length - 1).trim) }
          .filter { _.subtype.equalsIgnoreCase(config.subtype.get) }
          .map { _.sample }

        genePatientMatrix.retain { (pk, p) => samplesInfo.contains(pk.sample) }
      }

      genePatientMatrix.toMap
    }

    val all_samples = genePatientMatrix.map(_._1.sample).toSet
    /*
      println("===========================")
      all_samples.foreach { sample => println(sample) }
      println("===========================")
      */

    val geneList = new HashSet[Gene]
    network.genes.foreach(gene => {
      var count = 0
      for (sample <- all_samples) {
        if (genePatientMatrix contains (PolimorphismKey(gene, sample))) {
          count += 1
        }
      }

      if (count >= config.seedGenesMutations) {
        geneList.add(gene)
      }
    })

    (interactions, network, genePatientMatrix, all_samples, geneList.view.toSet)
  }

  def buildWalkers(geneList: Set[Gene], networkManager: MutualExclusivityNetworkManager) = {
    val walkers: Set[SubNetworkSelector] = geneList.map(gene =>
      //new Fungus(
      new Hi2iiFungus(
        startGene = gene,
        geneNumberVariable = 2 + StatUtils.getRandomPoisson(1),
        //geneNumberVariable = 3 + StatUtils.getRandomPoisson(0.5), // 0 - 60%, 1 - 30%, 2 - 10%
        //geneNumberVariable = 3,
        network = networkManager)).toSet

    walkers
  }

  def getCloseGenesRandomizedGenePatientMatrix(genePatientMatrix: Map[PolimorphismKey, Polimorphism], annotationMap: Map[Gene, AnnotationInfo]): Map[PolimorphismKey, Polimorphism] = {

    val toReturn = new HashMap[PolimorphismKey, Polimorphism]

    genePatientMatrix.keys.foreach { key =>

      val annotation = annotationMap.get(key.gene)
      if (annotation.isDefined) {
        val annotationInfo = annotation.get
        var newGene: Gene = annotationInfo.closeGenes(0).gene
        var ran = random.nextInt(annotationInfo.sumGenesSize)

        annotationInfo.closeGenes
          .takeWhile(_ => ran > 0)
          .foreach { annotation =>
            ran -= annotation.size
            newGene = annotation.gene
          }

        toReturn.put(PolimorphismKey(newGene, key.sample), Polimorphism(newGene.name))
      }
    }

    toReturn.view.toMap
  }

  def bootstraapSamples(alterationsBySample:Map[String, List[(PolimorphismKey, Polimorphism)]]): Map[PolimorphismKey, Polimorphism] = {
    
    val samples = alterationsBySample.keys.toArray
  
    List.tabulate(samples.size)( n => 
      alterationsBySample(samples(random.nextInt(samples.size))).map(e => (PolimorphismKey(e._1.gene, "Sample"+n),Polimorphism(e._2.geneSymbol)))
    ).flatten.toMap
  }
  
  def randomizeGeneNames(genePatientMatrix: Map[PolimorphismKey, Polimorphism], genesSet: Set[Gene]): Map[PolimorphismKey, Polimorphism] = {
    
    val originalGenes = genesSet.toList
    val randomizedGenes = random.shuffle(originalGenes)
    
    val old2new = originalGenes.zip(randomizedGenes).toMap
    
    genePatientMatrix
      .filter {e => genesSet contains e._1.gene}
      .map { e =>
        val newGene = old2new(e._1.gene)
        (PolimorphismKey(newGene, e._1.sample), Polimorphism(newGene.name))
      }.toMap
  }
  
  def randomizeBipartiteGenePatientMatrix(genePatientMatrix: Map[PolimorphismKey, Polimorphism], samplesSet: Set[String], genesSet: Set[Gene]): Map[PolimorphismKey, Polimorphism] = {

    print("Randomizing data... ")

    val start = System.currentTimeMillis()

    val samples = samplesSet.zipWithIndex.toMap
    val genes = genesSet.zipWithIndex.toMap

    val m = Array.fill(genes.size)(Array.fill(samples.size)(false))
    genes.par.foreach { g =>
      samples.foreach { s =>
        if (genePatientMatrix.contains(PolimorphismKey(g._1, s._1))) {
          m(g._2)(s._2) = true
        }
      }
    }
    println("Done creating m. [" + (System.currentTimeMillis() - start) + " ms]")

    val Np = 100 * genePatientMatrix.size
    val d = genePatientMatrix.size / (samples.size * genes.size)

    val N = Np * math.log((1 - d) * Np / 100) / (200 * (1 - d))

    val matrix = ArrayBuffer(genePatientMatrix.toList: _*)

    println("Done creating new Matrix. [" + (System.currentTimeMillis() - start) + " ms]")

    (0 to N.toInt).par.foreach { n =>

      val r1 = random.nextInt(matrix.size)
      val r2 = random.nextInt(matrix.size)

      val ab = matrix(r1)
      val cd = matrix(r2)

      if (r1 != r2 && ab._1.gene != cd._1.gene && ab._1.sample != cd._1.sample) {

        this.synchronized {
          val a = genes(ab._1.gene)
          val b = samples(ab._1.sample)
          val c = genes(cd._1.gene)
          val d = samples(cd._1.sample)
          if (!m(a)(d) && !m(c)(b)) {
            m(a)(d) = true
            m(c)(b) = true
            m(a)(b) = false
            m(c)(d) = false
            matrix(r1) = (PolimorphismKey(ab._1.gene, cd._1.sample), ab._2)
            matrix(r2) = (PolimorphismKey(cd._1.gene, ab._1.sample), cd._2)
          }
        }
      }
    }

    println("Done randomization. [" + (System.currentTimeMillis() - start) + " ms]")
    val toReturn = matrix.toMap
    println("Done. [" + (System.currentTimeMillis() - start) + " ms]")

    toReturn
  }

}

case class Annotation(val geneSymbol: String, val impact: String)
case class Polimorphism(val geneSymbol: String, val source: String = "any", val score: Double = 1, _type: String = "any")
case class PolimorphismKey(val gene: Gene, val sample: String)
case class SampleInfo(val sample: String, val subtype: String)

case class GeneAnnotation(gene: Gene, chrm: String, start: Int, end: Int) {
  lazy val size = end - start
}
case class AnnotationInfo(closeGenes: Seq[GeneAnnotation], sumGenesSize: Int)

case class Config(
  /*
   * Base parameters
   */
  iterations: Int,
  reinforcement: Double,
  forgetfulness: Double,
  refNetwork: Seq[String] = Seq("HT", "hiII14", "reactome"),
  useRank: Boolean = true,
  statistical: Boolean = false,
  minMutPerGene: Int = 3,
  seedGenesMutations: Int = 10,
  outputPrefix: String = "ME",
  input: Option[String] = None,
  outputGenes: Int = 100,
  subtype: Option[String] = None,
  subtypeFile: Option[File] = None,
  debug: Option[Seq[String]] = None,
  excludeGenes: Option[Seq[Gene]] = None,
  processors: Int = Runtime.getRuntime().availableProcessors() / 2,
  randomizeData: Boolean = false,

  /*
   * P-value calculation parameters
   */
  pvalueExperiments: Int = 1000,
  closestGenes: Int = 5,
  genomicdistnace: Int = 50000,
  pvaluetype:String = "distance",
  
  /*
   * Bootstraap calculation parameters
   */
  bootstraapExperiments:Int = 1000,
  minBootstraapSupport:Double = 0.95,
  positiveGeneSetLists: Seq[String] = Seq(),
  useCGC:Boolean=true,
  useNCG:Boolean=true,
  ppv:Double=0.6,
  
  /*
   * Print patern additional parameters
   */
  genes: List[String] = List())