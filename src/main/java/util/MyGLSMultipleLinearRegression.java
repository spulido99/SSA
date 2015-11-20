package util;

import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.regression.GLSMultipleLinearRegression;

public class MyGLSMultipleLinearRegression extends GLSMultipleLinearRegression {

	public double calculateResidualSumOfSquares() {
		final RealVector residuals = calculateResiduals();
        return residuals.dotProduct(residuals);
	}
	
}
