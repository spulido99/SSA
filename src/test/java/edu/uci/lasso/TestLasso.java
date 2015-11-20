package edu.uci.lasso;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

/**
 * This class demonstrates the usage of LassoFit 
 * 
 * @author Yasser Ganjisaffar (http://www.ics.uci.edu/~yganjisa/)
 */
public class TestLasso {

	public static void main(String[] args) throws Exception {
		
		/*
		 * Input data is in "src/test/resources/diabetes.data" file
		 * We initialize a reader to read this input file
		 */		
		BufferedReader reader = new BufferedReader(new InputStreamReader(TestLasso.class.getClassLoader().getResourceAsStream("diabetes.data")));
		
		/*
		 * The first line of the input file is the header which should be ignored.
		 * So, we read the first line
		 */
		String line = reader.readLine();
		
		/*
		 * Number of features (predictors) is determined based on the 
		 * number of columns in the header line
		 */
		String[] parts = line.split("\t");
		int featuresCount = parts.length - 1;
		
		/*
		 * Observations and targets are read and loaded from the input file
		 */
		List<float[]> observations = new ArrayList<float[]>();
		List<Float> targets = new ArrayList<Float>();
		while ((line = reader.readLine()) != null) {
			parts = line.split("\t");
			float[] curObservation = new float[featuresCount];
			for (int f = 0; f < featuresCount; f++) {
				curObservation[f] = Float.parseFloat(parts[f]);
			}
			observations.add(curObservation);
			targets.add(Float.parseFloat(parts[parts.length - 1]));
		}

		/*
		 * LassoFitGenerator is initialized
		 */
		LassoFitGenerator fitGenerator = new LassoFitGenerator();
		int numObservations = observations.size();
		fitGenerator.init(null, null);
		for (int i = 0; i < numObservations; i++) {
			fitGenerator.setObservationValues(i, observations.get(i));
			fitGenerator.setTarget(i, targets.get(i));
		}
		
		/*
		 * Generate the Lasso fit. The -1 arguments means that
		 * there would be no limit on the maximum number of 
		 * features per model
		 */
		LassoFit fit = fitGenerator.fit(-1);
		
		/*
		 * Print the generated fit
		 */
		System.out.println(fit);
	}
}