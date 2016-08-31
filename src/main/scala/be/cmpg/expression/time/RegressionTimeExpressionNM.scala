package be.cmpg.expression.time

import scala.collection.Set
import scala.collection.mutable.ArrayBuffer
import org.apache.commons.math3.linear.BlockRealMatrix
import org.apache.commons.math3.linear.SingularMatrixException
import org.apache.commons.math3.stat.correlation.Covariance
import be.cmpg.graph.Gene
import be.cmpg.graph.Interaction
import be.cmpg.graph.Network
import be.cmpg.graph.interaction.NodeCostNetworkManager
import util.MyGLSMultipleLinearRegression
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression
import be.cmpg.walk.SubNetworkSelector
import be.cmpg.utils.weightByFlatInitialProbability

class RegressionTimeExpressionNM(network: Network,
  timeExpressionMatrix: Map[String, Array[Array[Double]]], // gene -> Matrix[5 Samples][5 Time Points]
  toRegress: Array[Array[Double]],
  all_samples: Set[String], // all the samples names (here to not need to analyse all keys of the genePatientMatrix)
  mAS_perGene: Int, // minFreqGene: minimum frequency of mutation for a gene to be considered in the analysis
  pheromone: Double = 0.05,
  initialProb:Double = 0.5,
  evaporation: Double = 0.996,
  ranked: Boolean = false) extends NodeCostNetworkManager(network: Network, pheromone: Double, evaporation: Double, new weightByFlatInitialProbability(network,initialProb), ranked: Boolean) {

  override def scoreSubnetwork(subnetwork: Set[Interaction], selector: Option[SubNetworkSelector] = None): Double = {

    val genes = subnetwork.map(_.genes).flatten[Gene]

    if (genes.isEmpty)
      return 0.0

    val x = new ArrayBuffer[Array[Double]]
    val y = new ArrayBuffer[Double]

    for (i <- 0 until toRegress.size) {
      for (j <- 0 until toRegress(0).size) {

        val point = new ArrayBuffer[Double]

        genes.foreach(g => {
          val data = timeExpressionMatrix.get(g.name)
          if (data.isDefined) {
            val tmp = data.get(i)(j)
            point += tmp
          }
        })

        y += toRegress(i)(j)
        x += point.toArray
        //regression.addObservation(point.toArray, toRegress(i)(j))
      }
    }

    val regression = new OLSMultipleLinearRegression
    //val regression = new MyGLSMultipleLinearRegression
    //val regression = new LassoFitGenerator

    //regression.init(y.toArray, x.toArray);

    regression.newSampleData(y.toArray, x.toArray)
    /*
     * 
    val xM = new BlockRealMatrix(x.toArray)
    val covariance = new Covariance(xM.transpose()).getCovarianceMatrix().getData
    regression.newSampleData(y.toArray, x.toArray, covariance)
     */

    //val fit = regression.fit(-1);

    //val error = regression.calculateResidualSumOfSquares()
    val error = regression.calculateResidualSumOfSquares()

    if (error.isNaN() || error.isInfinite())
      return 0.0

    math.max(0.0, 1.0 - error)

  }

}