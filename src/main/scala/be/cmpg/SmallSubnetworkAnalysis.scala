package be.cmpg

import be.cmpg.cancer.MemoComparissonAnalysis
import be.cmpg.cancer.simulations.MutualExclusivitySimulations
import be.cmpg.cancer.PancancerAnalysis
import be.cmpg.cancer.MutualExclusivityAnalysis
import be.cmpg.cancer.MutualExclusivityPrintPattern
import be.cmpg.cancer.simulations.AnalizeSimulationResults
import be.cmpg.cancer.MutualExclusivityPrepareInput
import be.cmpg.cancer.PValueCalculator
import be.cmpg.cancer.BootstrapCalculator
import be.cmpg.cancer.pcawg.FunSeq2Analysis

object SmallSubnetworkAnalysis extends App {
  
  if (args.length == 0) {
    println("Usage: java -jar SSA.jar <Module> [options...]")
    println("Modules available: ")
    println("\t - Mutual Exclusivity: ME_input > ME (other: ME_input, ME_btstrp, ME_pattern)")
  } else {
    val module = args(0)
    val newArgs = args.drop(1)
    module match {
      case "ME" => MutualExclusivityAnalysis.main(newArgs)
      case "ME_input" => MutualExclusivityPrepareInput.main(newArgs)
      case "ME_pattern" => MutualExclusivityPrintPattern.main(newArgs)
      case "ME_pan" => PancancerAnalysis.main(newArgs)
      case "ME_funSeq2" => FunSeq2Analysis.main(newArgs)
      case "ME_memoComp" => MemoComparissonAnalysis.main(newArgs)
      case "ME_sim" => MutualExclusivitySimulations.main(newArgs)
      case "ME_sim_results" => AnalizeSimulationResults.main(newArgs)
      case "ME_pval" => PValueCalculator.main(newArgs)
      case "ME_btstrp" => BootstrapCalculator.main(newArgs)
    }
  }
  
}