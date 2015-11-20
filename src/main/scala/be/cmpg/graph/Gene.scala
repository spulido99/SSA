package be.cmpg.graph

import scala.collection.Map

case class Gene(val name: String, typ: String = "gene") {

  var data:Map[String, String] = _
  
  lazy val toStr = name+"["+typ+"]"
  
  override def equals(other: Any) = this.toStr.equals(other.toString())

  override def hashCode: Int = toStr.hashCode()
  
  override def toString = toStr
  
  def get(key:String) = data.get(key)
}