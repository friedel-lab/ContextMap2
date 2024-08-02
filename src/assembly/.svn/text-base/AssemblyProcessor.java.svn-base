package assembly;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Pattern;

public class AssemblyProcessor {

	public AssemblyProcessor() {
		
	}
	
	
	public void filterByLength(String assemblyFilePath, int minLength) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(assemblyFilePath)));
			String currentLine;
			Pattern spacePattern = Pattern.compile(" ");
			int currentLength;
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '>') {
					currentLength = Integer.valueOf(spacePattern.split(currentLine)[1]);
					if(currentLength >= minLength) {
						System.out.println(currentLine);
						System.out.println(br.readLine());
					}
				}
			}
			br.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public void getAssemblyStatistic(String assemblyFilePath) {
		try {
			int numberOfContigs;
			int numberOfContigsLarger100 = 0;
			int maxContigSize;
			long sumOfContigLengths = 0;
			double meanContigSize = 0;
			int n50Size = 0;
			int numberOfContigsLargerN50 = 0;
			int overallReads = 0;
			int coveredReads = 0;
			ArrayList<Integer> contigSizes = new ArrayList<Integer>();
			
			BufferedReader assemblyReader = new BufferedReader(new FileReader(new File(assemblyFilePath)));
			String currentLine;
			String[] splittedLine;
			int currentContigLength;
			while(assemblyReader.ready()) {
				currentLine = assemblyReader.readLine();
				if(currentLine.charAt(0) == '>') {
					splittedLine = currentLine.split(" ");
					currentContigLength = Integer.valueOf(splittedLine[1]);
					sumOfContigLengths += currentContigLength;
					if(currentContigLength > 100) numberOfContigsLarger100++;
					contigSizes.add(currentContigLength);
				}
			}
			Collections.sort(contigSizes);
			numberOfContigs = contigSizes.size();
			maxContigSize = contigSizes.get(contigSizes.size()-1);
			meanContigSize = (double)sumOfContigLengths/(double)numberOfContigs;
			long currentSumOfContigLengths = 0;
			for(int j = contigSizes.size() - 1; j >= 0; j --) {
				currentSumOfContigLengths += contigSizes.get(j);
				if(2 * currentSumOfContigLengths >= sumOfContigLengths) {
					n50Size = contigSizes.get(j);
					numberOfContigsLargerN50 = contigSizes.size() - j;
					break;
				}
			}
			System.out.println("Number of contigs: " + numberOfContigs);
			System.out.println("Mean size: " + meanContigSize);
			System.out.println("Max size: " + maxContigSize);
			System.out.println("N50 size: " + n50Size);
			System.out.println("Number of contigs >N50 size: " + numberOfContigsLargerN50);
			System.out.println("Sum of contig lenghts: " + sumOfContigLengths);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
}
