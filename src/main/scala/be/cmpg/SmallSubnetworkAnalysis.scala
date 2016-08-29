package be.cmpg

import be.cmpg.cancer.MemoAnalysis.MemoComparissonAnalysis
import be.cmpg.cancer.simulations.MutualExclusivitySimulations
import be.cmpg.cancer.pancancer.PancancerAnalysis
import be.cmpg.cancer.ME.MutualExclusivityAnalysis
import be.cmpg.cancer.functional.FunctionalMutualExclusivityAnalysis
import be.cmpg.cancer.simulations.AnalizeSimulationResults
import be.cmpg.cancer.ME.MutualExclusivityPrepareInput
import be.cmpg.statistical.PValueCalculator
import be.cmpg.statistical.BootstrapCalculator
import be.cmpg.cancer.functional.FunctionalMutualExclusivityPrepareInput
import be.cmpg.cluster.MushthofaAnalysis

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
      case "ME_pan" => PancancerAnalysis.main(newArgs)
      case "ME_memoComp" => MemoComparissonAnalysis.main(newArgs)
      case "ME_sim" => MutualExclusivitySimulations.main(newArgs)
      case "ME_sim_results" => AnalizeSimulationResults.main(newArgs)
      case "ME_pval" => PValueCalculator.main(newArgs)
      case "ME_btstrp" => BootstrapCalculator.main(newArgs)
      case "FME_input" => FunctionalMutualExclusivityPrepareInput.main(newArgs)
      case "FME" => FunctionalMutualExclusivityAnalysis.main(newArgs)
      case "P_val" => MushthofaAnalysis.main(newArgs)
    }
  }
  
}