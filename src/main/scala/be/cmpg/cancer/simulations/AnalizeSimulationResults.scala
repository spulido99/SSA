package be.cmpg.cancer.simulations

import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import scala.collection.JavaConversions._
import java.io.FileWriter
import be.cmpg.simulatedDataAnalysis.ROCPlotValue
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import scala.collection.mutable.HashMap
import java.io.File
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

object AnalizeSimulationResults extends App {
  
  var dir = "C:/Users/Spulido99/Documents/PhD_big/data/phac/"
  if (!(new File(dir)).exists())
    dir = ""
  
  
  /*
   * Parameter analysis
   */
  println("Parameter settings analysis...");

  if (false) {
    val filename = dir + "HT_hiII14simulation_output.parameters.txt"

    if (new File(filename).exists) {

      val writer = new FileWriter(filename.replaceAll(".txt", "") + ".analysis.txt")

      val parameterResults = new HashMap[Param, SummaryStatistics]

      
      val data = new CSVReader(new FileReader(filename), '\t')
          
      val byParam = new CSVReader(new FileReader(filename), '\t')
        .readAll()
        .drop(1)
        .groupBy { fields => Param(fields(1).toDouble, fields(2).toDouble) }
        .foreach { e =>
          val param = e._1

          val byRep = e._2.groupBy { fields => fields(0).toInt }

          val aucStat = new SummaryStatistics

          byRep.foreach { e =>
            val rocValues = e._2.map { f => new ROCPlotValue(f(4).toDouble, f(5).toDouble, f(6).toDouble, f(7).toDouble) }
            val auc = rocValues.foldLeft(0.0)((x, y) => x + y.sensitivity) / rocValues.size
            //println(param + " auc: " + auc)
            aucStat.addValue(auc)
          }

          writer.write(param.reinforcement + "\t" + param.forgetfulness + "\t" + aucStat.getMean + "\n")

        }
      writer.close()
    }
  }
  println("Done.")

  /*
   * The rest of the analysis
   */
  
  val fileNames = List(dir + "HT_hiII14simulation_output.subnetworksize.small.txt", dir + "HT_hiII14simulation_output.added.txt", dir + "HT_hiII14simulation_output.remove.txt")

  for (fileName <- fileNames) {
    val byN = new CSVReader(new FileReader(fileName), '\t')
      .readAll()
      .drop(1)
      .groupBy { fields => fields(2).toInt } // 0 NoiseLevel

    println(byN.size)
    
    var noiseLevels: List[String] = List()

    val writer = new FileWriter(fileName.replaceAll(".txt", "") + ".analysis.txt")
    byN.take(1).foreach(f => {
      noiseLevels = f._2.groupBy { _(0) }.keys.toList.sortBy { x => x }
      println(noiseLevels)
      writer.write("N")
      noiseLevels.foreach { x => 
        writer.write("\t" + x + "_phac_ppv\t" + x + "_phac_sensitivity\t" + x + "_phac_aspecificity") 
                  }
      writer.write("\n")
    })

    byN.keys.toList.sortBy { x => x }.foreach(keyN => {

      val byNoise = byN(keyN).groupBy { _(0) }

      writer.write(keyN.toString)
      noiseLevels.foreach(noise => {
        val data = byNoise.get(noise)

        if (data.isDefined) {
          
          
          val rocValues = data.get.map { f => new ROCPlotValue(f(3).toDouble, f(4).toDouble, f(5).toDouble, f(6).toDouble) }
          
          val statAspecificity = new DescriptiveStatistics
          val statSensitivity = new DescriptiveStatistics
          val statPPV = new DescriptiveStatistics
          
          rocValues.foreach { x => 
            statAspecificity.addValue(x.aspecificity) 
            statSensitivity.addValue(x.sensitivity)  
            statPPV.addValue(x.ppv)
          }
          
          writer.write("\t" + statPPV.getMean + "\t" + statSensitivity.getMean+ "\t" + statAspecificity.getMean)
              //+ "\t" + statAspecificity.getPercentile(25) + "\t" + statSensitivity.getPercentile(25)
              //+ "\t" + statAspecificity.getPercentile(75) + "\t" + statSensitivity.getPercentile(75))
        } else {
          writer.write("\t" + 1.0 + "\t" + 1.0 + "\t" + 0.0)
        }
      })
      writer.write("\n")
    })
    writer.close()
  }

}