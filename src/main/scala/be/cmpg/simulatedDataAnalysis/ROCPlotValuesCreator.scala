package be.cmpg.simulatedDataAnalysis

import java.io.FileWriter

class ROCPlotValuesCreator(selected: Set[String], positives: Set[String], all: Set[String]) {

  val value = {

    val negatives = all.--(positives)
    
    val truePositives = positives.&(selected)
    val trueNegatives = (negatives).&(all.--(selected))
    val falsePositives = selected.&(negatives)
    val falseNegatives = positives.--(selected)

    new ROCPlotValue(truePositives.size.toDouble, falsePositives.size.toDouble, trueNegatives.size.toDouble, falseNegatives.size.toDouble)
    
    //val sensitivity = truePositives.size.toDouble / (truePositives.size.toDouble + falseNegatives.size.toDouble)
    //val positivePredictiveValue = truePositives.size.toDouble / (truePositives.size.toDouble + falsePositives.size.toDouble)

    //(positivePredictiveValue , sensitivity)
  }
  
  def get = value

}
