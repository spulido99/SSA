package be.cmpg

import be.cmpg.cancer.MemoComparissonAnalysis
import be.cmpg.cancer.simulations.MutualExclusivitySimulations
import be.cmpg.cancer.PancancerAnalysis
import be.cmpg.cancer.MutualExclusivityAnalysis
import be.cmpg.cancer.MutualExclusivityPrintPattern
import be.cmpg.cancer.simulations.AnalizeSimulationResults
import be.cmpg.cancer.MutualExclusivityPrepareInput
import be.cmpg.cancer.PValueCalculator

object SmallSubnetworkAnalysis extends App {  
  
  if (args.length == 0) {
    println("Usage: java -jar SSA.jar <Module> [options...]")
    println("Modules available: ")
    println("\t - Mutual Exclusivity: ME_input > ME (other: ME_pattern, ME_pan, ME_memoComp, ME_sim, ME_pval)")
  } else {
    val module = args(0)
    val newArgs = args.drop(1)
    module match {
      case "ME" => MutualExclusivityAnalysis.main(newArgs)
      case "ME_input" => MutualExclusivityPrepareInput.main(newArgs)
      case "ME_pattern" => MutualExclusivityPrintPattern.main(newArgs)
      case "ME_pan" => PancancerAnalysis.main(newArgs)
      case "ME_memoComp" => MemoComparissonAnalysis.main(newArgs)
      case "ME_sim" => MutualExclusivitySimulations.main(newArgs)
      case "ME_sim_results" => AnalizeSimulationResults.main(newArgs)
      case "ME_pval" => PValueCalculator.main(newArgs)
    }
  }
  
}