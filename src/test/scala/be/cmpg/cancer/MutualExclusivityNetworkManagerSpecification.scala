package be.cmpg.cancer

import org.junit.runner.RunWith
import org.specs2.runner.JUnitRunner
import org.specs2.mutable.Specification
import be.cmpg.graph.Network
import be.cmpg.graph.Interaction
import be.cmpg.graph.Gene

@RunWith(classOf[JUnitRunner])
class MutualExclusivityNetworkManagerSpecification extends Specification {

  "The Mutual Exclusivity Network Manager" should {
    "give score of this matrix should be the expected score" in {
      val network = new Network(interactions = Set())

      /*
     * Score:
     * 
     * Mutation matrix
     *      g1   g2   g3   g4  
     * S_1   x                       
     * S_2   x         x    x 
     * S_3   x                
     * S_4   x                
     * S_5   x    x          
     * S_6        x    x       
     * S_7        x          
     * S_8        x           
     * S_9             x      
     * S_10                 x
     *  .
     *  .
     *  .     
     * S_30                           // samples without mutation in those genes
     */

      val genePatientMatrix: Map[PolimorphismKey, Polimorphism] = Map(
        PolimorphismKey(Gene("g1"), "S_1",0) -> Polimorphism("g1"),
        PolimorphismKey(Gene("g1"), "S_2",0) -> Polimorphism("g1"),
        PolimorphismKey(Gene("g3"), "S_2",0) -> Polimorphism("g3"),
        PolimorphismKey(Gene("g4"), "S_2",0) -> Polimorphism("g4"),
        PolimorphismKey(Gene("g1"), "S_3",0) -> Polimorphism("g1"),
        PolimorphismKey(Gene("g1"), "S_4",0) -> Polimorphism("g1"),
        PolimorphismKey(Gene("g1"), "S_5",0) -> Polimorphism("g1"),
        PolimorphismKey(Gene("g2"), "S_5",0) -> Polimorphism("g2"),
        PolimorphismKey(Gene("g2"), "S_6",0) -> Polimorphism("g2"),
        PolimorphismKey(Gene("g3"), "S_6",0) -> Polimorphism("g3"),
        PolimorphismKey(Gene("g2"), "S_7",0) -> Polimorphism("g2"),
        PolimorphismKey(Gene("g2"), "S_8",0) -> Polimorphism("g2"),
        PolimorphismKey(Gene("g3"), "S_9",0) -> Polimorphism("g3"),
        PolimorphismKey(Gene("g4"), "S_10",0) -> Polimorphism("g4")) ++ (11 to 30).map {id => PolimorphismKey(Gene("gOther"), "S_"+id,0) -> Polimorphism("gOther")}.toMap

      /* Scores:
       *            S1    S2  S3  S4   S5  S6  S7  S8  S9  S10
     * s(g1) = sqrt( 1 + 1/3 + 1 + 1 +1/2)                    = sqrt(3.8333)  = 1.957
     * s(g2) = sqrt(                      1/2 + 1 + 1)        = sqrt(2.5)     = 2.5
     * s(g3) = sqrt(1) = 1.73
     * s(g4) = sqrt(1) = 1.73
     * 
     * S_2 y S_5 only contribute once to the score of g1. 
     * S_2 is mutated in 3 genes, so the score increases for g1 is only 1/3
     * S_5 is mutated in 2 genes, so the score increases for g1 is only 1/2
     * 
     * For g3, S_2 and S_6 does not add to the score as it have already given a score to g1 and g2 respectively.
     * For g4, S_2 does not add to the score as it have already given to g1. 
     *      
     */

      val networkManager = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix,
        minimumSamplesAltered = 0,
        evaporation = 0.996)

      val subnetwork = Set(
        Interaction(Gene("g1"), Gene("g2")),
        Interaction(Gene("g2"), Gene("g3")),
        Interaction(Gene("g3"), Gene("g4")))

      val expectedScore = (math.sqrt(1.0 + 1.0/3.0 + 1.0 + 1.0 + 1.0/2.0) + math.sqrt(1.0 + 1.0 + 1.0/2.0) + math.sqrt(1.0) + math.sqrt(1.0)) /4
      
      networkManager.scoreSubnetwork(subnetwork) must beEqualTo(expectedScore)
    }
    /*
    "better distributed networks should have better scores" in {
      val network = new Network(interactions = Set())

      
     * Score:
     * 
     * Mutation matrix
     *      g1   g2   g3   g4       g1   g2   g3   g4  
     * S_1   x                       x                       
     * S_2   x                       x                      
     * S_3   x                       x
     * S_4   x                            x
     * S_5   x                   <        x
     * S_6        x                       x
     * S_7        x                            x 
     * S_8        x                            x
     * S_9             x                       x
     * S_10            x                           x
     * S_11                x                       x
     * S_12                x                       x
     

      val genePatientMatrix1: Map[(String, String), Polimorphism] = Map(
        ("g1", "S_1") -> Polimorphism("g1"),
        ("g1", "S_2") -> Polimorphism("g1"),
        ("g1", "S_3") -> Polimorphism("g1"),
        ("g1", "S_4") -> Polimorphism("g1"),
        ("g1", "S_5") -> Polimorphism("g1"),
        ("g2", "S_6") -> Polimorphism("g2"),
        ("g2", "S_7") -> Polimorphism("g2"),
        ("g2", "S_8") -> Polimorphism("g2"),
        ("g3", "S_9") -> Polimorphism("g3"),
        ("g3", "S_10") -> Polimorphism("g3"),
        ("g4", "S_11") -> Polimorphism("g4"),
        ("g4", "S_12") -> Polimorphism("g4"))

      val genePatientMatrix2: Map[(String, String), Polimorphism] = Map(
        ("g1", "S_1") -> Polimorphism("g1"),
        ("g1", "S_2") -> Polimorphism("g1"),
        ("g1", "S_3") -> Polimorphism("g1"),
        ("g2", "S_4") -> Polimorphism("g2"),
        ("g2", "S_5") -> Polimorphism("g2"),
        ("g2", "S_5") -> Polimorphism("g2"),
        ("g2", "S_6") -> Polimorphism("g2"),
        ("g3", "S_7") -> Polimorphism("g3"),
        ("g3", "S_8") -> Polimorphism("g3"),
        ("g3", "S_9") -> Polimorphism("g3"),
        ("g4", "S_10") -> Polimorphism("g4"),
        ("g4", "S_11") -> Polimorphism("g4"),
        ("g4", "S_12") -> Polimorphism("g4"))

      val networkManager1 = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix1,
        minimumSamplesAltered = 0)

      val networkManager2 = new MutualExclusivityNetworkManager(
        network = network,
        genePatientMatrix = genePatientMatrix2,
        minimumSamplesAltered = 0)

      val subnetwork = Set(
        Interaction(Gene("g1"), Gene("g2")),
        Interaction(Gene("g2"), Gene("g3")),
        Interaction(Gene("g3"), Gene("g4")))

      println("BadSubNetwork: "+ networkManager1.scoreSubnetwork(subnetwork))
      println("GoodSubNetwork: "+ networkManager2.scoreSubnetwork(subnetwork))
      
      networkManager1.scoreSubnetwork(subnetwork) must beLessThan(networkManager2.scoreSubnetwork(subnetwork) )
    }*/
  }

  /*
  "The Mutual Exclusivity Network Manager" should {

    val network = new Network(
      interactions = Set(
        Interaction(Gene("exclusive1"), Gene("exclusive2"), probability = 1),
        Interaction(Gene("exclusive1"), Gene("notmutated1"), probability = 1),
        Interaction(Gene("exclusive1"), Gene("notexclusive1"), probability = 1),
        Interaction(Gene("notmutated1"), Gene("notexclusive2"), probability = 1),
        Interaction(Gene("notexclusive1"), Gene("notexclusive"), probability = 1),
        Interaction(Gene("exclusive2"), Gene("exclusive3"), probability = 1),
        Interaction(Gene("exclusive2"), Gene("exclusive4"), probability = 1)))

    val genePatientMatrix: Map[(String, String), String] = Map(
      // p1 have only exclusive1
      (("exclusive1", "p1.1") -> "M"),
      (("exclusive1", "p1.2") -> "M"),
      (("exclusive1", "p1.3") -> "M"),
      // p2 have only exclusive2
      (("exclusive2", "p2.1") -> "M"),
      (("exclusive2", "p2.2") -> "M"),
      // p3 have only exclusive3
      (("exclusive3", "p3.1") -> "M"),
      (("exclusive3", "p3.2") -> "M"),
      // p4 have only exclusive4
      (("exclusive4", "p4.1") -> "M"),
      (("exclusive4", "p4.2") -> "M"),
      // notexclusive 1&2 have mutations for pacients 1,2 & 3
      (("notexclusive1", "p1.1") -> "M"),
      (("notexclusive1", "p1.2") -> "M"),
      (("notexclusive2", "p1.1") -> "M")
      )

    val all_samples = genePatientMatrix.keys.map(_._2).toSet

    val networkManager = new MutualExclusivityNetworkManager(
      network = network,
      genePatientMatrix = genePatientMatrix,
      mAS_perGene = 1,
      all_samples = all_samples,
      evaporation = 0.996)

    "give zero score to an empty subnetwork" in {
      val subnetwork: Set[Interaction] = Set()
      networkManager.scoreSubnetwork(subnetwork) must beEqualTo(0)
    }

    "give zero score to a size one subnetwork" in {
      val subnetwork: Set[Interaction] = Set(Interaction(Gene("exclusive1"), Gene("exclusive1")))
      networkManager.scoreSubnetwork(subnetwork) must beEqualTo(0)
    }

    "give zero score to a subnetwork with non mutated genes" in { // FIXME should it be like this??? What if all subnetowrks have genes without mutations???
      val subnetwork: Set[Interaction] = Set(Interaction(Gene("exclusive1"), Gene("notmutated1")))
      networkManager.scoreSubnetwork(subnetwork) must beEqualTo(0)
    }

    "give a positive score to a 2 gene subnetwork" in {
      val subnetwork: Set[Interaction] = Set(Interaction(Gene("exclusive1"), Gene("exclusive2")))
      networkManager.scoreSubnetwork(subnetwork) must be_>(0.0)
    }

    "give a higher positive score to bigger exclusive subnetworks" in {
      val subnetwork_s2: Set[Interaction] = Set(Interaction(Gene("exclusive1"), Gene("exclusive2")))
      val subnetwork_s3: Set[Interaction] = Set(
        Interaction(Gene("exclusive1"), Gene("exclusive2")),
        Interaction(Gene("exclusive2"), Gene("exclusive3")))

      val score_2 = networkManager.scoreSubnetwork(subnetwork_s2)
      val score_3 = networkManager.scoreSubnetwork(subnetwork_s3)

      score_3 must be_>(score_2)
    }

    "give a higher score to subnetworks with more mutual exclusive samples" in {
      val subnetwork_with5patients: Set[Interaction] = Set(Interaction(Gene("exclusive1"), Gene("exclusive2")))
      val subnetwork_with4patients: Set[Interaction] = Set(Interaction(Gene("exclusive3"), Gene("exclusive2")))

      val score_with5patients = networkManager.scoreSubnetwork(subnetwork_with5patients)
      val score_with4patients = networkManager.scoreSubnetwork(subnetwork_with4patients)

      score_with5patients must be_>(score_with4patients)
    }

	    "give a higher score to subnetworks without non-exclusive genes" in {
      val subnetwork_withNonExclusive: Set[Interaction] = Set(
          Interaction(Gene("exclusive1"), Gene("exclusive2")),
          Interaction(Gene("exclusive1"), Gene("notexclusive1")))
      val subnetwork_withOnlyExclusive: Set[Interaction] = Set(
          Interaction(Gene("exclusive1"), Gene("exclusive2")))

      val score_withNonExclusive = networkManager.scoreSubnetwork(subnetwork_withNonExclusive)
      val score_withOnlyExclusive = networkManager.scoreSubnetwork(subnetwork_withOnlyExclusive)

      score_withOnlyExclusive must be_>(score_withNonExclusive)
    }
  }
  *
  */
}