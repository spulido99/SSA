package be.cmpg.walk.neighbourhoodScoring

class MutationScoresGenerator(MutationList:List[String]) {

  val genesWithMutations = MutationList.map(entry => entry.split("\t")(0)).toSet
  
// this method calculates the scores based on the quotient of non_synonymous mutations and all mutations
  def getScoresQuotient:scala.collection.mutable.HashMap[String,String] ={
  val scoreMap = new scala.collection.mutable.HashMap[String,String]()
  genesWithMutations.foreach(gene =>{
    var synonymousMutations = 0
    var nonSynonymousMutations =0
    MutationList.foreach(mutation =>{
      val splitted = mutation.split("\t")
      if(splitted(0).equals(gene)){
        if(splitted(1).equals("No"))
          synonymousMutations = synonymousMutations+1
          else nonSynonymousMutations = nonSynonymousMutations+1}
    })
    if(!(synonymousMutations+nonSynonymousMutations == 0)){
  scoreMap.put(gene,(nonSynonymousMutations.toDouble/(synonymousMutations.toDouble+nonSynonymousMutations.toDouble)).toString +"\t"+  nonSynonymousMutations.toString +"\t"+  synonymousMutations.toString)}
  
  else
  scoreMap.put(gene,"0")
  })
  return scoreMap
  }
  
  // this method calculates the scores by simply counting the mutations
  def getScoresCount:scala.collection.mutable.HashMap[String,String] ={
  val scoreMap = new scala.collection.mutable.HashMap[String,String]()
  genesWithMutations.foreach(gene =>{
    var synonymousMutations = 0
    var nonSynonymousMutations =0
    MutationList.foreach(mutation =>{
      val splitted = mutation.split("\t")
      if(splitted(0).equals(gene)){
        if(splitted(1).equals("No"))
          synonymousMutations = synonymousMutations+1
          else nonSynonymousMutations = nonSynonymousMutations+1}
    })
    if(!(synonymousMutations+nonSynonymousMutations == 0)){
  scoreMap.put(gene,(nonSynonymousMutations.toDouble+synonymousMutations.toDouble).toString+"\t"+ nonSynonymousMutations.toString +"\t"+ synonymousMutations.toString)}
  
  else
   scoreMap.put(gene,"0")
  })
  return scoreMap
  }
  
  // this method calculates the scores by deviding the number of non-synonymous mutations by the number of mutations
  def getScoresEvolution:scala.collection.mutable.HashMap[String,String] ={
  val scoreMap = new scala.collection.mutable.HashMap[String,String]()
  genesWithMutations.foreach(gene =>{
    var synonymousMutations = 0
    var nonSynonymousMutations =0
    MutationList.foreach(mutation =>{
      val splitted = mutation.split("\t")
      if(splitted(0).equals(gene)){
        if(splitted(1).equals("No"))
          synonymousMutations = synonymousMutations+1
          else nonSynonymousMutations = nonSynonymousMutations+1}
    })
    if(!(synonymousMutations+nonSynonymousMutations == 0)){
  scoreMap.put(gene,(nonSynonymousMutations.toDouble/synonymousMutations.toDouble).toString +"\t"+ nonSynonymousMutations.toString +"\t"+ synonymousMutations.toString)}
  
  else
   scoreMap.put(gene,"0")
  })
  return scoreMap
  }
  
}