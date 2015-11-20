package resultAnalysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class RocAnalysis {

	private static class AnalysisInfo {
		double phac_aspecificity;
		double phac_sensitivity;
		double neighbourhood_aspecificity;
		double neighbourhood_sensitivity;
		double stat_aspecificity;
		double stat_sensitivity;
	}
	
	public static void main(String[] args) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader("human10_output.txt"));
		String line;
		
		AnalysisInfo[] info = new AnalysisInfo[16027];
		for (int i = 0; i < info.length; i ++) {
			info[i] = new AnalysisInfo();
		}
		
		int n = 0;
		boolean first = true;
		while ((line = br.readLine()) != null) {
			if (first) {
				first = false;
			} else {
				String[] data = line.split("\t");
				int pos = 0;
				n = Integer.parseInt(data[pos++]);
				int i = Integer.parseInt(data[pos++]);
				info[i].phac_aspecificity += Double.parseDouble(data[pos++]);
				info[i].phac_sensitivity += Double.parseDouble(data[pos++]);
				info[i].neighbourhood_aspecificity += Double.parseDouble(data[pos++]);
				info[i].neighbourhood_sensitivity += Double.parseDouble(data[pos++]);
				info[i].stat_aspecificity += Double.parseDouble(data[pos++]);
				info[i].stat_sensitivity += Double.parseDouble(data[pos++]);
			}
		}
		br.close();
		
		BufferedWriter writer = new BufferedWriter(new FileWriter("human10_output_average.txt"));
		
		for (int i = 0; i < info.length; i ++) {
			StringBuilder sb = new StringBuilder(i);
			sb.append("\t").append(info[i].phac_aspecificity / n);
			sb.append("\t").append(info[i].phac_sensitivity / n);
			sb.append("\t").append(info[i].neighbourhood_aspecificity / n);
			sb.append("\t").append(info[i].neighbourhood_sensitivity / n);
			sb.append("\t").append(info[i].stat_aspecificity / n);
			sb.append("\t").append(info[i].stat_sensitivity / n);
			sb.append("\n");
			
			writer.append(sb);
		}
		writer.flush();
		writer.close();
		
	}
	
}
