package be.cmpg.simulatedDataAnalysis

import math.sqrt

class ROCPlotValue(val _TP:Double, val _FP:Double, val _TN:Double, val _FN:Double) {
  
  lazy val sensitivity = _TP / (_TP + _FN)
  
  lazy val specificity = _TN / (_FP + _TN)
  
  lazy val aspecificity = 1 - specificity
  
  lazy val ppv = _TP / (_TP + _FP)
  
  /*
   * https://en.wikipedia.org/wiki/Matthews_correlation_coefficient
   */
  lazy val mcc = (_TP*_TN - _FP*_FN)/sqrt((_TP+_FP)*(_TP+_FN)*(_TN+_FP)*(_TN+_FN))
}