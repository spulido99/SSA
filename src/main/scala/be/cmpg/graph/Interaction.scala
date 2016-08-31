package be.cmpg.graph

case class Interaction(from: Gene, to: Gene,
  typ: String = "any",
  direction: String = "undirected",
  regulatory: String = "regulatory",
  // Probability of interaction was not being used so disabled for now as to not confuse.
//  probability: Double = 0d,
  evidence:Set[String] = Set()) {

  lazy val genes = Set(from, to)
  
  // It does not make sense to compare the names of the nodes in order to decide on the direction of the interaction
  
/*  lazy val toStr = {
    if (direction != "undirected" || from.name < to.name) {
      from.name + "\t" + typ + "\t" + to.name
    } else {
      // always return the same str if undirected
      to.name + "\t" + typ + "\t" + from.name
    }
  }*/
  
  lazy val toStr ={
    from.name + "\t" + typ + "\t" + to.name +"\t" + direction
  }
  
  def contains(gene:Gene){from==Gene || to ==Gene}
  
  def listGenes = List(from,to)
  
  override def equals(other: Any) = this.toStr.equals(other.toString())

  override def hashCode: Int = toStr.hashCode()

  override def toString = toStr

}