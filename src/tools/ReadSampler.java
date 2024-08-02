package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.regex.Pattern;

public class ReadSampler {

	
	public ReadSampler() {
		
	}
	
	
	public void sampleReads(String readFilePath, int sampleSize, String outputFilePath) {
		try {
			Date date = new Date();
			System.out.println(String.format("[%s]\tGetting number of reads.",date.toLocaleString()));
			int numberOfReads = getNumberOfReads(readFilePath);
			
			if(sampleSize > numberOfReads) {
				System.err.println("sample size > #reads in input file. Returning...");
				System.exit(1);
			}
	
			
			date = new Date();
			System.out.println(String.format("[%s]\tGenerating sample.",date.toLocaleString()));
			HashSet<Integer> samples = generateSample(sampleSize,numberOfReads);
			
			
			
			date = new Date();
			System.out.println(String.format("[%s]\tWriting reads.",date.toLocaleString()));
			int readPosition = 0;
			BufferedReader br  = new BufferedReader(new FileReader(new File(readFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			String currentLine;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) == '>') {
					if(samples.contains(readPosition)) {
						pw.println(currentLine);
						pw.println(br.readLine());
					}
					readPosition++;
				}
			}
			br.close();
			pw.close();
			System.out.println(String.format("[%s]\tDone.",date.toLocaleString()));
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public int getNumberOfReads(String readFilePath) throws Exception {
		int numberOfReads = 0;
		BufferedReader br = new BufferedReader(new FileReader(new File(readFilePath)));
		String currentLine;
		while((currentLine = br.readLine()) != null) {
			if(currentLine.charAt(0) == '>') numberOfReads++;
		}
		br.close();
		return numberOfReads;
	}
	
	public HashSet<Integer> generateSample(int sampleSize, int numberOfReads) throws Exception {
		HashSet<Integer> samples = new HashSet<Integer>();
		Random random = new Random();
		
		while(samples.size() < sampleSize) {
			samples.add(random.nextInt(numberOfReads));
		}
		return samples;
		
	}
	
	
	public void getSampleStatistic(String samplingFolderPath, double minCov) {
		try {
			HashMap<String,Species> genbankid2species = new HashMap<String,Species>();
			ArrayList<Integer> samples = new ArrayList<Integer>();
			File[] sampleFolders = new File(samplingFolderPath).listFiles();
			BufferedReader contaminationReader;
			BufferedReader confidenceReader;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			
			for(File sampleFolder : sampleFolders) {
				int currentSample = Integer.valueOf(sampleFolder.getName());
				samples.add(currentSample);
				
				//read top10 contamination file
				contaminationReader = new BufferedReader(new FileReader(new File(sampleFolder.getAbsolutePath() + "/contamination.txt.top10")));
				while((currentLine = contaminationReader.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					if(!genbankid2species.containsKey(splittedLine[0])) {
						genbankid2species.put(splittedLine[0], new Species(splittedLine[1],Integer.valueOf(splittedLine[4])));
					}
					genbankid2species.get(splittedLine[0]).addSample(currentSample, Integer.valueOf(splittedLine[2]));
				}
				contaminationReader.close();
				
				//read confidence file
				confidenceReader = new BufferedReader(new FileReader(new File(sampleFolder.getAbsolutePath() + "/confidence.txt")));
				confidenceReader.readLine();
				while((currentLine = confidenceReader.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					if(genbankid2species.containsKey(splittedLine[0])) {
						genbankid2species.get(splittedLine[0]).addConfidence(currentSample,Double.valueOf(splittedLine[4]));
					}
				}
				confidenceReader.close();
			}
			
			//filter by coverage
			Species species;
			int currentGenomeSize;
			int currentReadCounts;
			double currentCoverage;
			boolean minCovReached;
			
			HashMap<Integer,Integer> sample2readCounts;
			ArrayList<String> toDelete = new ArrayList<String>();
			if(minCov > 0) {
				for(String speciesName : genbankid2species.keySet()) {
					species = genbankid2species.get(speciesName);
					currentGenomeSize = species.getGenomeSize();
					sample2readCounts = species.getSample2readcounts();
					minCovReached = false;
					for(int sample : sample2readCounts.keySet()) {
						currentReadCounts = sample2readCounts.get(sample);
						currentCoverage = (double)currentReadCounts/(double)currentGenomeSize;
						species.addCoverage(sample, currentCoverage);
						if(currentCoverage >= minCov) {
							minCovReached = true;
						}
					}
					
					if(!minCovReached)
						toDelete.add(speciesName);
				}
				
				for(String speciesName : toDelete)
					genbankid2species.remove(speciesName);
			}
			
			
			Collections.sort(samples);
			System.out.print("species_name\ttrivial_name\tgenome_size");
			for(int sample : samples) {
				System.out.print("\t" + sample + "_r\t" + sample + "_conf\t" + sample + "_cov");
			}
			System.out.println();
			Species currentSpecies;
			
			HashMap<Integer,Double> sample2confidence;
			HashMap<Integer,Double> sample2coverage;
			for(String speciesName : genbankid2species.keySet()) {
				currentSpecies = genbankid2species.get(speciesName);
				System.out.print(speciesName + "\t" + currentSpecies.getTrivialName() + "\t" + currentSpecies.getGenomeSize());
				sample2readCounts = currentSpecies.getSample2readcounts();
				sample2confidence = currentSpecies.getSample2confidence();
				sample2coverage = currentSpecies.getSample2coverage();
				for(int sample : samples) {
					if(sample2readCounts.containsKey(sample)) {
						System.out.print("\t" + sample2readCounts.get(sample) + "\t" + sample2confidence.get(sample) + "\t" + sample2coverage.get(sample));
					}
					else
						System.out.print("\t0\t0\t0");
				}
				System.out.println();
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private class Species {
		
		private String trivialName;
		private int genomeSize;
		private HashMap<Integer,Integer> sample2readcounts;
		private HashMap<Integer,Double> sample2confidence;
		private HashMap<Integer,Double> sample2coverage;
		
		public Species(String trivialName,int genomeSize) {
			this.trivialName = trivialName;
			this.genomeSize = genomeSize;
			this.sample2readcounts = new HashMap<Integer,Integer>();
			this.sample2confidence = new HashMap<Integer,Double>();
			this.sample2coverage = new HashMap<Integer,Double>();
		}
		
		public void addSample(int sampleSize,int readCount) {
			this.sample2readcounts.put(sampleSize, readCount);
		}
		
		public String getTrivialName() {
			return this.trivialName;
		}
		
		public int getGenomeSize() {
			return this.genomeSize;
		}
		
		public HashMap<Integer,Double> getSample2confidence() {
			return this.sample2confidence;
		}
		
		public void addConfidence(int sample, double confidence) {
			this.sample2confidence.put(sample,confidence);
		}
		
		public HashMap<Integer,Double> getSample2coverage() {
			return this.sample2coverage;
		}
		
		public void addCoverage(int sample, double coverage) {
			this.sample2coverage.put(sample,coverage);
		}
		
		public HashMap<Integer,Integer> getSample2readcounts() {
			return this.sample2readcounts;
		}
	}
}
