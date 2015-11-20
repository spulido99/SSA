package be.cmpg.utils

import scala.util.Random

object StatUtils {

  def getRandomPoisson(mean: Double) = {
    val L = math.exp(-mean)
    var k = 0
    var p = 1.0

    do {
      k += 1
      p = p * Random.nextDouble
    } while (p > L)

    k - 1
  }

}