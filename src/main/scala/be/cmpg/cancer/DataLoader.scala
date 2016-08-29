package be.cmpg.cancer

import java.io.File
import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.util.Try
import scala.collection.JavaConversions.asScalaBuffer
import java.util.LinkedList
import scala.io.Source
import scala.util.Random
import java.io.InputStreamReader
import java.io.FileInputStream
import be.cmpg.graph.Gene
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import be.cmpg.graph.Network

object DataLoader {

  // Set up random generator
  val random = new Random(System.nanoTime())

  // Load true positive genes for assessing the PPV of your network. For this CGC and NCG are being used
  def cgcGenes(inputFolder: String) = {
    new CSVReader(new InputStreamReader(new FileInputStream(inputFolder + "CancerGeneCensus20150718.csv"))) // CGC: Cancer Gene Census
      .readAll()
      .map(_(0))
      .map(Gene(_))
      .toSet
  }

  def ncgGenes(inputFolder: String) = {

    new CSVReader(new InputStreamReader(new FileInputStream(inputFolder + "NCG5.0.txt")), '\t') // NCG 5.0: Network of Cancer Genes
      .readAll()
      .map(_(0))
      .map(Gene(_))
      .toSet

  }

  // This is a set of genes which is sometimes used to filter CNV's. When used, genes with CNV's not in this list get filtered out
  def zackGenes(inputFolder: String): Set[String] = {
    Source.fromInputStream(new FileInputStream(inputFolder + "ZackEtAl_CNAGenes.txt")).getLines.toSet
  }

  def loadExpression(cnvdata: Map[String, File], inputFolder: String, acceptedCorrelationQVal: Double = 0.05, gisticThreshold: Int = 1, maxQtyCopyNumber: Int = 500, GenesConsidered: List[String] = List(), filterBasedOnKnownGenes: Boolean = false): Map[PolymorphismKey, Polymorphism] = {

    val gistic = cnvdata.get("gistic")
    val peaks = cnvdata.get("cnv_peaks")
    val expression = cnvdata.get("exp")
    val correlation = cnvdata.get("corr")

    // Gistic output is needed as it calls CNV's and gives thresholds for each CNV
    if (gistic.isEmpty)
      throw new RuntimeException("If using CNV gistic folder should exist.")

    val expressionR: Map[String, Double] =
      // If correlation files are readily available, use them
      if (correlation.isDefined) {
        println("CNV: loading correlation")

        new CSVReader(new FileReader(correlation.get), '\t')
          .readAll()
          .filter { fields => Try { fields(3).toDouble }.getOrElse(1.0) <= acceptedCorrelationQVal } // Column (field) 3 is the qvalue
          .map { fields =>
            val gene = fields(0)
            val corr = Try { fields(1).toDouble }.getOrElse(0.0)
            (gene, corr)
          }.toMap

        // If correlation files are not readily available, calculate correlation from expression data and peak data
      } else if (expression.isDefined && peaks.isDefined) {

        println("CNV: Calculating correlation")

        val allexp = readFromGeneSampleMatrix(expression.get)
        val allcnvs = readFromGeneSampleMatrix(peaks.get)

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

          // calculate correlation,r, of CNV and expression data for gene g
          val r = ((X zip Y) map { case (x, y) => (x - uX) * (y - uY) }).sum / (sX * sY)

          if (!r.isNaN() && !r.isInfinite())
            expressionR.add((g, r))
        }

        val max = expressionR.map(_._2).max

        // use only positive correlations and scale so the correlation value is 1     
        expressionR.toList.filter(_._2 > 0.0).map { case (g, r) => (g, r / max) }.toMap

      } else {
        throw new RuntimeException("No Correlation or files to calculate correlation available.")
      }

    println("CNV: Reading peaks")

    val nonsamples = List("Gene Symbol", "Locus ID", "Cytoband")

    // Load amplifications and deletions from the gistic files which have at least 99% confidence
    val gisticAmpGenes = loadGisticConfidenceGenes(gistic.get + "/amp_genes.conf_99.txt")

    val gisticDelGenes = loadGisticConfidenceGenes(gistic.get + "/del_genes.conf_99.txt")

    // Either use all GISTIC confident amplifications and deletions or only those which are allready known cancer genes
    val confGenes = if (filterBasedOnKnownGenes) { (gisticAmpGenes ++ gisticDelGenes).toSet intersect (zackGenes(inputFolder) ++ ncgGenes(inputFolder) ++ cgcGenes(inputFolder)) }
    else { (gisticAmpGenes ++ gisticDelGenes).toSet }

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

    // Take only the confident genes defined previously which also correlate with expression data and which have a gistic threshold of at least the defined value
    val all = gisticOutput.map { fields =>
      var count = fields.size - samples.size - 1
      val geneName = fields(0)
      samples.map { sample =>
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
            // Location is not yet implemented here !!
            if (GenesConsidered.isEmpty) { Some(PolymorphismKey(Gene(geneName), sample, 0), Polymorphism(geneName, if (value > 0) "amp" else "del", geneExpressionCorr.get)) }
            else {
              if (GenesConsidered.contains(geneName)) { Some(PolymorphismKey(Gene(geneName), sample, 0), Polymorphism(geneName, if (value > 0) "amp" else "del", geneExpressionCorr.get)) }
              else {
                None
              }
            }
          }
        }
      }
    }.flatten.flatten.toList

    // Identify the hypermutator samples
    val hypermutators = all.groupBy(_._1.sample).filter(_._2.size > maxQtyCopyNumber).keys.toSet

    // Filter out all CNV's belonging to hypermutators
    all.filter { e => !hypermutators.contains(e._1.sample) }.toMap
  }

  def loadGisticConfidenceGenes(fileName: String) = {

    Source.fromFile(fileName, "latin1")
      .getLines
      .drop(4)
      .map { x => x.split("\t") }
      .flatten
      .filter { x => x.length() > 0 }
  }


  // Reads a maf file and puts it in a map. Maf can contain functional scores.
  def loadMaf(mafFile: File, dropLines: Int = 2, geneNameCol: Int = 0, sampleCol: Int = 15, locationCol: Int = 5, substringSample: Boolean = true, mutationTypeCol: Int = 8, toExclude: List[String] = List("Silent", "Intron", "LOH"), maxQtyMutations: Int = 500, GenesConsidered: List[String] = List(),functionalData:Boolean=false): Map[PolymorphismKey, Polymorphism] = {

    val all = Source.fromFile(mafFile, "latin1")
      .getLines
      .drop(dropLines)
      .filter(line => !line.startsWith("\t"))
      .map { _.split("\t") }
      .map { fields =>

        val geneName = fields(geneNameCol)
        var sample = fields(sampleCol)
        val location = fields(locationCol).toInt
        val Score = if (functionalData){fields.last.toDouble}
        else{1d}
        if (substringSample)
          sample = sample.substring(0, 15)

        val _type = fields(mutationTypeCol)

        if (toExclude.contains(_type)) {
          None
        } else {
          if (GenesConsidered.isEmpty) {
            Some(PolymorphismKey(Gene(geneName), sample.trim.replaceAll("\\.", "-"), location), Polymorphism(geneName, "mut", Score, _type))
          } else {
            if (GenesConsidered.contains(geneName)) { Some(PolymorphismKey(Gene(geneName), sample.trim.replaceAll("\\.", "-"), location), Polymorphism(geneName, "mut", Score, _type)) }
            else {
              None
            }
          }
        }
      }
      .flatten
      .toList

    // Hypermutator patients are left out
    val hypermutators = all.groupBy(_._1.sample).filter(_._2.size > maxQtyMutations).keys.toSet

    all.filter { e => !hypermutators.contains(e._1.sample) }.toMap

  }

  // Load mutations and CNV's from TCGA2012 files for MEMo comparison
  def loadFromTcga2012Files(inputFolder: String) = {

    println("Reading TCGA 2012: Somatic Mutations maf file")
    val toReturnM = loadMaf(new File("src/test/resources/tcga/BRCA2012/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf"))

    println("Reading TCGA 2012: Calculating CNV correlation with expression")

    val toReturnE = loadExpression(Map(
      "cnv_peaks" -> new File("src/test/resources/tcga/BRCA2012/TCGA_BRCApub_cnv_peaks.txt"),
      "exp" -> new File("src/test/resources/tcga/BRCA2012/BRCA.exp.547.med.txt"),
      "cnv_peaks" -> new File("src/test/resources/tcga/BRCA2012/TCGA_BRCApub_cnv_peaks_th.txt")), inputFolder)

    println("Reading TCGA 2012: Done.")
    toReturnM ++ toReturnE
  }

  // Read CNV and expression data files in which there is data for CNV/expression per gene per patient
  def readFromGeneSampleMatrix(file: File, dropFromHeader: Int = 1, dropFromFields: Int = 0, valueToIgnore: String = "null") = {
    // Hashmap of (sample,gene),value needs to be returned
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
            if (v != valueToIgnore && !v.isEmpty() && v != "NA")
              toReturn.put((geneName, sample.trim.replaceAll("\\.", "-")), v.toDouble)
            count += 1
          })
        }
      })

    toReturn.view.toMap
  }

  // Load mutationMatrix file which may or may not contain functional data as well
  def loadMutationMatrixFiles(input: String, functionalData: Boolean = false): HashMap[PolymorphismKey, Polymorphism] = {

    val genePatientMatrix = new HashMap[PolymorphismKey, Polymorphism]
    new CSVReader(new FileReader(input), '\t').readAll().foreach { fields =>
      val sample = fields(0)
      fields.drop(1).foreach { gene =>
        // Location not yet implemented here !
        val key = if (functionalData) { PolymorphismKey(Gene(gene.split(",")(0)), sample, 0) }
        else { PolymorphismKey(Gene(gene), sample, 0) }
        val value = if (functionalData) { Polymorphism(gene.split(",")(0), score = gene.split(",")(1).toDouble) }
        else { Polymorphism(gene) }

        genePatientMatrix.put(key, value)
      }
    }

    genePatientMatrix
  }
  
  // Loads the needed data for analysis.
    def loadData(config: Config, FunctionalData: Boolean) = {
    val excludedGenesList = if (config.excludedGenes==""){List[String]()} 
    else{Source.fromFile(new File(config.excludedGenes)).getLines.toList}
    val interactions = NetworkFactory.loadNetwork(config.refNetwork, config.inputFolder,excludedGenesList)
    val network = new Network(interactions)
    val networkFolder = config.inputFolder

    println("Interactions: " + interactions.size)

    val genePatientMatrix = {

      println("Loading mutation matrix file...")
      val genePatientMatrix = if (FunctionalData) {
        this.loadMutationMatrixFiles(config.input.get,true)
      } else {
        this.loadMutationMatrixFiles(config.input.get)
      }
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
        // Location not yet implemented here !
        if (genePatientMatrix contains (PolymorphismKey(gene, sample, 0))) {
          count += 1
        }
      }

      if (count >= config.seedGenesMutations) {
        geneList.add(gene)
      }
    })

    (interactions, network, genePatientMatrix, all_samples, geneList.view.toSet)
  }
}