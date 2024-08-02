package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Locale;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.regex.Pattern;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;

public class Statistic {

	private static HashSet<String> vertebrateChr;
	private static HashSet<String> vertebrateRDNA;
	private static  HashSet<String> microbialChr;
	private static  HashSet<String> virusChr;
	private static  HashSet<String> yeastChr;
	
	
	public Statistic() {
		//init hashs
		vertebrateChr = new HashSet<String>();
		microbialChr = new HashSet<String>();
		microbialChr.add("microbes_0");
		for(int i = 1; i <= 22; i++) {
			vertebrateChr.add(String.format("chr%s",i));
			microbialChr.add(String.format("microbes_%s",i));
		}
		vertebrateChr.add("chrX");
		vertebrateChr.add("chrY");
		vertebrateChr.add("chrM");
		
		vertebrateRDNA = new HashSet<String>();
		vertebrateRDNA.add("chrRDNA");
		
		virusChr = new HashSet<String>();
		virusChr.add("virus_0");
		
		yeastChr = new HashSet<String>();
		yeastChr.add("chrI");
		yeastChr.add("chrII");
		yeastChr.add("chrIII");
		yeastChr.add("chrIV");
		yeastChr.add("chrV");
		yeastChr.add("chrVI");
		yeastChr.add("chrVII");
		yeastChr.add("chrVIII");
		yeastChr.add("chrIX");
		yeastChr.add("chrYeastX");
		yeastChr.add("chrXI");
		yeastChr.add("chrXII");
		yeastChr.add("chrXIII");
		yeastChr.add("chrXIV");
		yeastChr.add("chrXV");
		yeastChr.add("chrXVI");
		yeastChr.add("chrYeastM");
	}
	
	public static void main(String args[]) {
		
		//init hashs
		vertebrateChr = new HashSet<String>();
		microbialChr = new HashSet<String>();
		microbialChr.add("microbes_0");
		for(int i = 1; i <= 22; i++) {
			vertebrateChr.add(String.format("chr%s",i));
			microbialChr.add(String.format("microbes_%s",i));
		}
		vertebrateChr.add("chrX");
		vertebrateChr.add("chrY");
		vertebrateChr.add("chrM");
		
		vertebrateRDNA = new HashSet<String>();
		vertebrateRDNA.add("chrRDNA");
		
		virusChr = new HashSet<String>();
		virusChr.add("virus_0");
		
		yeastChr = new HashSet<String>();
		yeastChr.add("chrI");
		yeastChr.add("chrII");
		yeastChr.add("chrIII");
		yeastChr.add("chrIV");
		yeastChr.add("chrV");
		yeastChr.add("chrVI");
		yeastChr.add("chrVII");
		yeastChr.add("chrVIII");
		yeastChr.add("chrIX");
		yeastChr.add("chrYeastX");
		yeastChr.add("chrXI");
		yeastChr.add("chrXII");
		yeastChr.add("chrXIII");
		yeastChr.add("chrXIV");
		yeastChr.add("chrXV");
		yeastChr.add("chrXVI");
		yeastChr.add("chrYeastM");
		
		if(args[0].equals("getErrorRate")) {
			String samFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
			}
			getErrorRate(samFilePath);
		}
		
		//needs the default bowtie output in the input
		if(args[0].equals("getErrorRateFromBowtieOutput")) {
			String bwtFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-bwt")) {
					bwtFilePath = args[++i];
				}
			}
			getErrorRateFromBowtieOutput(bwtFilePath);
		}
		
		if(args[0].equals("getReadIdsFromBowtieOutput")) {
			String bwtFilePath = null;
			int maxMismatches = -1;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-bwt")) {
					bwtFilePath = args[++i];
				}
				if(args[i].equals("-mm")) {
					maxMismatches = Integer.valueOf(args[++i]);
				}
			}
			getReadIdsFromBowtieOutput(bwtFilePath,maxMismatches);
		}
		
		if(args[0].equals("getLengthDistribution")) {
			String fastaFilePath = null;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					fastaFilePath = args[++i];
				}
			}
			getLengthDistribution(fastaFilePath);
		}
		
		if(args[0].equals("getErrorRateTable")) {
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			double minCoverage = Double.MIN_VALUE;
			boolean useMDflag = false;
			HashSet<String> referenceChromosomes = new HashSet<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-microbialcontent")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("-reference")) {
					String[] references = args[++i].split(",");
					for(String ref : references)
						referenceChromosomes.add(ref);
					continue;
				}
				if(args[i].equals("-mincov")) {
					minCoverage = Double.valueOf(args[++i]);
				}
				if(args[i].equals("--mdflag")) {
					useMDflag = true;
				}
				
			}
			generateErrorRateTableFromConcatenatedGenomes(samFilePath,microbialContentFilePath,genomeIndexDirPath,useMDflag,referenceChromosomes,null);
		}
		
		if(args[0].equals("averageErrorRateTable")) {
			ArrayList<String> filePaths = new ArrayList<String>();
			int maxMismatches = -1;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					String[] splittedPaths = args[++i].split(",");
					filePaths.addAll(Arrays.asList(splittedPaths));
				}
				if(args[i].equals("-mismatches")) {
					maxMismatches = Integer.valueOf(args[++i]);
				}
			}
			averageErrorRateTables(filePaths,maxMismatches);
		}
		
		if(args[0].equals("getDistanceTable")) {
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean mergeYeastChr = true;
			HashSet<String> speciesFilterNames = new HashSet<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-microbialcontent")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("-filter")) {
					String[] names = args[++i].split(",");
					speciesFilterNames.addAll(Arrays.asList(names));
				}
				if(args[i].equals("--noyeastmerge")) {
					mergeYeastChr = false;
				}
			}
			generateDistanceTableFromConcatenatedGenomes(samFilePath, microbialContentFilePath,genomeIndexDirPath,speciesFilterNames,mergeYeastChr);
		}
		
		if(args[0].equals("getDistanceTableFromUnconcatenatedGenomes")) {
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean mergeYeastChr = true;
			HashSet<String> speciesFilterNames = new HashSet<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-microbialcontent")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("-filter")) {
					String[] names = args[++i].split(",");
					speciesFilterNames.addAll(Arrays.asList(names));
				}
				if(args[i].equals("--noyeastmerge")) {
					mergeYeastChr = false;
				}
			}
			generateDistanceTableFromUnconcatenatedGenomes(samFilePath, microbialContentFilePath,genomeIndexDirPath,speciesFilterNames,mergeYeastChr);
		}
		
		if(args[0].equals("filterDistanceTable")) {
			String tableFilePath = null;
			HashSet<String> speciesNames = new HashSet<String>();
			double maxDistance = Integer.MAX_VALUE;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					tableFilePath = args[++i];
				}
				if(args[i].equals("-species")) {
					String[] splittedNames = args[++i].split(",");
					for(int j = 0; j < splittedNames.length; j++) {
						speciesNames.add(splittedNames[j]);
					}
				}
				if(args[i].equals("-maxdistance")) {
					maxDistance = Double.valueOf(args[++i]);
				}
				
			}
			filterDistanceTable(tableFilePath, speciesNames, maxDistance);
		}
		
		if(args[0].equals("getConfidenceScores")) {
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean mergeYeastChr = true;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-microbialcontent")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("--noyeastmerge")) {
					mergeYeastChr = false;
				}
			}
			generateConfidenceScoresFromConcatenatedGenomes(samFilePath, microbialContentFilePath,genomeIndexDirPath,mergeYeastChr,null);
		}
		
		if(args[0].equals("getConfidenceScoresFromUnconcatenatedGenomes")) {
			String samFilePath = null;
			String microbialContentFilePath = null;
			String genomeIndexDirPath = null;
			boolean mergeYeastChr = true;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("-microbialcontent")) {
					microbialContentFilePath = args[++i];
				}
				if(args[i].equals("-idx")) {
					genomeIndexDirPath = args[++i];
				}
				if(args[i].equals("--noyeastmerge")) {
					mergeYeastChr = false;
				}
			}
			generateConfidenceScoresFromUnconcatenatedGenomes(samFilePath, microbialContentFilePath,genomeIndexDirPath,mergeYeastChr);
		}
		
		if(args[0].equals("averageConfidenceScores")) {
			ArrayList<String> filePaths = new ArrayList<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					String[] splittedPaths = args[++i].split(",");
					filePaths.addAll(Arrays.asList(splittedPaths));
				}
			}
			averageConfidenceScores(filePaths);
		}
		
		
		
		if(args[0].equals("readDistribution")) {
			String samFilePath = null;
			Boolean splitMicrobes = false;
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
				}
				if(args[i].equals("--splitMicrobes")) {
					splitMicrobes = true;
				}
			}
			getMappingStats(samFilePath,splitMicrobes);
		}
		
		if(args[0].equals("compareDistributions")) {
			String referenceDistPath = null;
			ArrayList<String> otherDistPaths = new ArrayList<String>();
			ArrayList<String> names = new ArrayList<String>();
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-ref")) {
					referenceDistPath = args[++i];
					continue;
				}
				if(args[i].equals("-compare")) {
					String[] paths = args[++i].split(",");
					otherDistPaths.addAll(Arrays.asList(paths));
					continue;
				}
				if(args[i].equals("-names")) {
					names.addAll(Arrays.asList(args[++i].split(",")));
					continue;
				}
				
			}
			compareDist(referenceDistPath,otherDistPaths,names);
		}
		
		//compares local with global read scores. Therefore the final mapping file as well as the file
		//'all_resolved_local_contexts.txt.sorted.best.matchings' are needed.
		if(args[0].equals("compareContextScores")) {
			String globalMappingFilePath = null;
			String localMappingFilePath = null;
			String outputDirPath = null;
			String indexFolderPath = null;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-globalmapping")) {
					globalMappingFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-localmapping")) {
					localMappingFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx")) {
					indexFolderPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputDirPath = args[++i];
					continue;
				}
				
			}
			compareContextScores(globalMappingFilePath,localMappingFilePath,indexFolderPath,outputDirPath);
		}
		
		if(args[0].equals("getStartPositionsFromUnconcatenatedGenomes")) {
			String samFilePath = null;
			String indexDirPath = null;
			String outputFilePath = null;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx")) {
					indexDirPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
			}
			getAbsoluteStartPositionsFromUnconcatenatedGenomes(samFilePath,indexDirPath,outputFilePath);
		}
		
		if(args[0].equals("getStartPositionsFromConcatenatedGenomes")) {
			String samFilePath = null;
			String indexDirPath = null;
			String outputFilePath = null;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-sam")) {
					samFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-idx")) {
					indexDirPath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
			}
			getAbsoluteStartPositionsFromConcatenatedGenomes(samFilePath,indexDirPath,outputFilePath);
		}
		
		if(args[0].equals("getAccuracyFromTable")) {
			String tableFilePath = null;
			String indexDirPath = null;
			String outputFilePath = null;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-table")) {
					tableFilePath = args[++i];
					continue;
				}
				
			}
			extractAccuracyFromTable(tableFilePath);
		}
		
		if(args[0].equals("getWorkload")) {
			String mpstatFilePath = null;
			String outputFilePath = null;
			int maxcores = Integer.MIN_VALUE;
			
			for(int i = 1; i < args.length; i++) {
				if(args[i].equals("-i")) {
					mpstatFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-o")) {
					outputFilePath = args[++i];
					continue;
				}
				if(args[i].equals("-maxcores")) {
					maxcores = Integer.valueOf(args[++i]);
					continue;
				}
				
			}
			buildWorkloadStatistic(mpstatFilePath, outputFilePath, maxcores);
		}
		
		
		
	}
	
	
	private static void getAbsoluteStartPositionsFromUnconcatenatedGenomes(String samFilePath, String indexDirPath, String outputPath) {
		try {
			HashMap<String,Microbe> idx2microbe = new HashMap<String,Microbe>();
			File[] indexFiles = new File(indexDirPath).listFiles();
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			String indexId;
			String shortId;
			BufferedReader br;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						indexId = splittedLine[0];
						splittedId = barPattern.split(splittedLine[1]);
						
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
						Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[2]),0);
						idx2microbe.put(indexId,tmpMicrobe);
					}
					br.close();
					
				}
			}
			
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			String readId;
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				indexId = splittedLine[2];
				if(idx2microbe.containsKey(indexId)) {
					pw.println(String.format("%s::%s\t%s",splittedLine[0],idx2microbe.get(indexId).getId(),splittedLine[3]));
				}
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void getAbsoluteStartPositionsFromConcatenatedGenomes(String samFilePath, String indexDirPath, String outputPath) {
		try {
			//now generate interval trees on all contained genomes
			BufferedReader br;
			HashMap<String,IntervalTree<Microbe>> index2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
			File[] indexFiles = new File(indexDirPath).listFiles();
			String currentLine;
			 String[] splittedLine;
			 String[] splittedId;
			 String shortId;
			 Pattern tabPattern = Pattern.compile("\t");
			 Pattern barPattern = Pattern.compile("\\|");
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					IntervalTree<Microbe> genomeTree = new IntervalTree<Microbe>();
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						splittedId = barPattern.split(splittedLine[0]);
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
						
						Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
						genomeTree.add(tmpMicrobe);
					}
					br.close();
					index2genomeTree.put(indexFile.getName().substring(0, indexFile.getName().lastIndexOf('.')), genomeTree);
				}
			}
			
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputPath)));
			IntervalTree<Microbe> currentGenomeTree;
			Microbe currentMicrobe;
			int relativeStart;
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				if(index2genomeTree.containsKey(splittedLine[2])) {
					currentGenomeTree = index2genomeTree.get(splittedLine[2]);
					relativeStart = Integer.valueOf(splittedLine[3]);
					currentMicrobe = currentGenomeTree.getIntervalsSpanning(relativeStart, new ArrayList<Microbe>()).get(0);
					pw.println(String.format("%s::%s\t%s",splittedLine[0],currentMicrobe.getId(),relativeStart - currentMicrobe.getStart()));
				}
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * in case we have mappings to microbes, virus or human microbiome we have local contexts overlapping several species!
	 * therefore we have to check every multi mapping for this and add it to the global context multi hits if necessary
	 */
	
	private static void compareContextScores(String globalMappingFilePath, String localMappingFilePath, String indexFolderPath, String outputDirPath) {
		try {
			
			BufferedReader br;
			String currentLine;
			String[] splittedLine;
			String[] splittedId;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern barPattern = Pattern.compile("\\|");
			String shortId;
			
			//parse global mapping file
			BufferedReader globalReader = new BufferedReader(new FileReader(new File(globalMappingFilePath)));
			PrintWriter globalWriter = new PrintWriter(new FileWriter(new File(outputDirPath + "/global_scores.txt")));
			
			Pattern doublePointPattern = Pattern.compile(":");
			double bestHitScore;
			double secondBestHitScore;
			while((currentLine = globalReader.readLine()) != null) {
				splittedLine =  tabPattern.split(currentLine);
				//found multi mapping
				if(splittedLine.length >= 12) {
					bestHitScore = Double.MIN_VALUE;
					secondBestHitScore = Double.MIN_VALUE;
					for(int i = 12; i < splittedLine.length; i++) {
						if(splittedLine[i].contains("S1:f:")) {
							bestHitScore = Double.valueOf(doublePointPattern.split(splittedLine[i])[2]);
							continue;
						}
						if(splittedLine[i].contains("S2:f:")) {
							secondBestHitScore = Double.valueOf(doublePointPattern.split(splittedLine[i])[2]);
						}
					}
					
					if(bestHitScore != Double.MIN_VALUE && secondBestHitScore != Double.MIN_VALUE) {
						globalWriter.println(bestHitScore - secondBestHitScore);
						if(bestHitScore - secondBestHitScore < 0) {
							System.err.println("WARNING: Found score < 0: " + currentLine);
						}
					}
				}
			}
			globalReader.close();
		
			if(localMappingFilePath != null) {
				//generate interval trees on all contained genomes
				HashMap<String,IntervalTree<Microbe>> index2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
				File[] indexFiles = new File(indexFolderPath).listFiles();
				for(File indexFile : indexFiles) {
					if(indexFile.getName().contains(".idx")) {
						br = new BufferedReader(new FileReader(indexFile));
						IntervalTree<Microbe> genomeTree = new IntervalTree<Microbe>();
						while(br.ready()) {
							currentLine = br.readLine();
							splittedLine = tabPattern.split(currentLine);
							splittedId = barPattern.split(splittedLine[0]);
							
							if(splittedId.length < 3)
								shortId = splittedId[0];
							else
								shortId = splittedId[3];
								
							Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
							genomeTree.add(tmpMicrobe);
						}
						br.close();
						index2genomeTree.put(indexFile.getName().substring(0, indexFile.getName().lastIndexOf('.')), genomeTree);
					}
				}
				
				
				//parse local mapping file
				BufferedReader localReader = new BufferedReader(new FileReader(new File(localMappingFilePath)));
				PrintWriter localWriter = new PrintWriter(new FileWriter(new File(outputDirPath + "/local_scores.txt")));
				String bestHitChr;
				String secondBestHitChr;
				
				int bestHitStart;
				int secondBestHitStart;
				
				Microbe bestHitMicrobe;
				Microbe seconbestHitMicrobe;
				
				IntervalTree<Microbe> currentTree;
				boolean overlapsGenomes;
				while((currentLine = localReader.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					//found multi mapping
					if(splittedLine.length >= 12) {
						overlapsGenomes = false;
						bestHitChr = splittedLine[2];
						secondBestHitChr = splittedLine[13];
						
						if(bestHitChr.equals(secondBestHitChr) && index2genomeTree.containsKey(bestHitChr)) {
							bestHitStart = Integer.valueOf(splittedLine[4]);
							secondBestHitStart = Integer.valueOf(splittedLine[14]);
							currentTree =  index2genomeTree.get(bestHitChr);
							
							bestHitMicrobe = currentTree.getIntervalsSpanning(bestHitStart, new ArrayList<Microbe>()).get(0);
							seconbestHitMicrobe = currentTree.getIntervalsSpanning(secondBestHitStart, new ArrayList<Microbe>()).get(0);
							
							if(!bestHitMicrobe.equals(secondBestHitChr)) {
								overlapsGenomes = true;
							}
								
						}
						
						if(!overlapsGenomes)
							localWriter.println(splittedLine[11]);
						
						//TODO write to third file here
						else
							globalWriter.println(splittedLine[11]);
					}
				}
				
				localReader.close();
				localWriter.close();
			}
			globalWriter.close();
		}
		catch(Exception e) { 
			e.printStackTrace();
		}
	}
	
	
	
	
	
	private static void getMappingStats(String samFilePath,Boolean splitMicrobes) {
		try {
			int vertebrateReads = 0;
			int vertebrateRdnaReads = 0;
			int microbialReads = 0;
			int virusReads = 0;
			int yeastReads = 0;
			int unknownReads = 0;
			int microbialReads_Index1 = 0;
			int microbialReads_Index2 = 0;
			
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			String currentLine;
			String[] splittedLine;
			String chr;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern underscorePattern = Pattern.compile("_");
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				
				if(vertebrateChr.contains(chr)) {
					vertebrateReads++;
					continue;
				}
				if(vertebrateRDNA.contains(chr)) {
					vertebrateRdnaReads++;
					continue;
				}
				if(microbialChr.contains(chr)) {
					microbialReads++;
					if(Integer.valueOf(underscorePattern.split(chr)[1]) <= 11)
						microbialReads_Index1++;
					else
						microbialReads_Index2++;
					
					continue;
				}
				if(virusChr.contains(chr)) {
					virusReads++;
					continue;
				}
				if(yeastChr.contains(chr)) {
					yeastReads++;
					continue;
				}
				unknownReads++;
			}
			br.close();
			
			System.out.println("type\tread count");
			System.out.println(String.format("vertebrate\t%s",vertebrateReads));
			System.out.println(String.format("rDNA\t%s",vertebrateRdnaReads));
			if(!splitMicrobes)
				System.out.println(String.format("microbe\t%s",microbialReads));
			else {
				System.out.println(String.format("microbe_1\t%s",microbialReads_Index1));
				System.out.println(String.format("microbe_2\t%s",microbialReads_Index2));
			}
			System.out.println(String.format("virus\t%s",virusReads));
			System.out.println(String.format("yeast\t%s",yeastReads));
			System.out.println(String.format("unknown\t%s",unknownReads));
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void compareDist(String refDistPath, ArrayList<String> compareDistPaths,ArrayList<String> names) {
		try {
			//read ref dist
			BufferedReader br = new BufferedReader(new FileReader(new File(refDistPath)));
			HashMap<String,Integer> refDist = new HashMap<String,Integer>();
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.length() == 0)
					continue;
				splittedLine = tabPattern.split(currentLine);
				refDist.put(splittedLine[0], Integer.valueOf(splittedLine[1]));
			}
			br.close();
			
			String currentName;
			HashMap<String,HashMap<String,Double>> names2relativeDist = new HashMap<String,HashMap<String,Double>>();
			HashMap<String,Double> currentRelativeDist;
			for(int i = 0; i < compareDistPaths.size(); i++) {
				currentName = names.get(i);
				currentRelativeDist = new HashMap<String,Double>();
				br = new BufferedReader(new FileReader(new File(compareDistPaths.get(i))));
				String species;
				int refCount;
				int currentCount;
				while(br.ready()) {
					currentLine = br.readLine();
					if(currentLine.length() == 0)
						continue;
					splittedLine = tabPattern.split(currentLine);
					species = splittedLine[0];
					currentCount = Integer.valueOf(splittedLine[1]);
					if(refDist.containsKey(species)) {
						currentRelativeDist.put(species, (double)((double)currentCount/(double)refDist.get(species)));
					}
				}
				br.close();
				names2relativeDist.put(currentName, currentRelativeDist);
			}
			
			System.out.print("species\tref");
			for(String name : names)
				System.out.print("\t" + name);
			System.out.print("\n");
			
			for(String species : refDist.keySet()) {
				System.out.print(String.format("%s\t1.0",species));
				for(String name : names) {
					currentRelativeDist = names2relativeDist.get(name);
					if(currentRelativeDist.containsKey(species))
						System.out.print(String.format("\t%s",currentRelativeDist.get(species)));
					else
						System.out.print("\t0.0");
				}
				System.out.print("\n");
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void getErrorRate(String samFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(samFilePath)));
			HashMap<Integer,Integer> missmatchCounts = new HashMap<Integer,Integer>();
			int maxMissmatches = -1;
			String currentLine;
			String[] splittedLine;
			StringTokenizer st;
			String currentToken = null;
			int currentMissmatchCount;
			while(br.ready()) {
				currentLine = br.readLine();
				st = new StringTokenizer(currentLine,"\t");
				currentMissmatchCount = -1;
				while(st.hasMoreTokens()) {
					currentToken = st.nextToken();
					if(currentToken.contains("NM:i:")) {
						currentMissmatchCount = Integer.valueOf(currentToken.substring(currentToken.lastIndexOf(":") + 1));
					}
				}
				
				if(currentMissmatchCount == -1) {
					currentMissmatchCount = Integer.valueOf(currentToken);
				}
				
				if(missmatchCounts.containsKey(currentMissmatchCount))
					missmatchCounts.put(currentMissmatchCount, missmatchCounts.get(currentMissmatchCount) + 1);
				else
					missmatchCounts.put(currentMissmatchCount,1);
				if(currentMissmatchCount > maxMissmatches)
					maxMissmatches = currentMissmatchCount;
			}
			br.close();
			
			System.out.println("missmatches\tread count");
			for(int i = 0; i <= maxMissmatches; i++) {
				if(missmatchCounts.containsKey(i))
					System.out.println(String.format("%s\t%s",i,missmatchCounts.get(i)));
				else
					System.out.println(String.format("%s\t%s",i,0));
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void getErrorRateFromBowtieOutput(String bwtFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(bwtFilePath)));
			HashMap<Integer,Integer> missmatchCounts = new HashMap<Integer,Integer>();
			int maxMissmatches = -1;
			String currentLine;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedLine;
			StringTokenizer st;
			String currentToken = null;
			int currentMissmatchCount;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				currentMissmatchCount = 0;
				if(splittedLine.length >= 8) {
					currentMissmatchCount = commaPattern.split(splittedLine[7]).length;
				}
				
				if(missmatchCounts.containsKey(currentMissmatchCount))
					missmatchCounts.put(currentMissmatchCount, missmatchCounts.get(currentMissmatchCount) + 1);
				else
					missmatchCounts.put(currentMissmatchCount,1);
				if(currentMissmatchCount > maxMissmatches)
					maxMissmatches = currentMissmatchCount;
			}
			br.close();
			
			System.out.println("missmatches\tread count");
			for(int i = 0; i <= maxMissmatches; i++) {
				if(missmatchCounts.containsKey(i))
					System.out.println(String.format("%s\t%s",i,missmatchCounts.get(i)));
				else
					System.out.println(String.format("%s\t%s",i,0));
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void getReadIdsFromBowtieOutput(String bwtFilePath, int maxMismatches) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(bwtFilePath)));
			String currentLine;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedLine;
			int currentMismatchCount;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				currentMismatchCount = 0;
				if(splittedLine.length >= 8) {
					currentMismatchCount = commaPattern.split(splittedLine[7]).length;
				}
				
				if(currentMismatchCount <= maxMismatches)
					System.out.println(splittedLine[0]);
					
			}
			br.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private static void getLengthDistribution(String fastaFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fastaFilePath)));
			TreeMap<Integer,MutableInt> lengthDistribution = new TreeMap<Integer,MutableInt>();
			String currentLine;
			int overallReadCount = 0;
			int overallLengths = 0;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) != '>') {
					overallReadCount++;
					overallLengths += currentLine.length();
					if(lengthDistribution.containsKey(currentLine.length()))
						lengthDistribution.get(currentLine.length()).increment();
					else
						lengthDistribution.put(currentLine.length(),new MutableInt(1));
				}
			}
			br.close();
			System.out.println("length\tcount");
			
			double averageLength = (double)overallLengths/(double)overallReadCount;
			double deviation = 0;
			for(int length : lengthDistribution.keySet()) {
				System.out.println(length + "\t" + lengthDistribution.get(length));
				deviation += (double)(lengthDistribution.get(length).intValue() * ((length - averageLength) * (length - averageLength)))/(double)overallReadCount;
			}
			System.out.println();
			
			System.out.println("average length: " + averageLength);
			System.out.println("sd: " + Math.sqrt(deviation));
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public static HashMap<String,HashMap<Integer,Double>> generateErrorRateTableFromConcatenatedGenomes(String samFilePath, String microbialContentFilePath,String genomeIndexDirPath, boolean useMDflag, HashSet<String> referenceChromosomes,HashMap<String,String> id2trivialName) {
		try {
			HashMap<String,HashMap<Integer,Double>> species2errorRate = new HashMap<String,HashMap<Integer,Double>>();
			HashMap<String,String> species2trivialName = new HashMap<String,String>();
			BufferedReader br;
			String currentLine;
			String[] splittedLine;
			String[] splittedId;
			String shortId;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern barPattern = Pattern.compile("\\|");
			String tmpChr;
			
			if(id2trivialName != null) {
				for(String id : id2trivialName.keySet()) {
					species2errorRate.put(id, new HashMap<Integer,Double>());
					species2trivialName.put(id, id2trivialName.get(id));
				}
			}
			
			else if(microbialContentFilePath != null) {
				br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
				//skipping the header
				br.readLine();
				while(br.ready()) {
					currentLine = br.readLine();
					splittedLine = tabPattern.split(currentLine);
					species2errorRate.put(splittedLine[0],new HashMap<Integer,Double>());
					species2trivialName.put(splittedLine[0],splittedLine[1]);
					
				}
				br.close();
			}
			
			//also add an entry for the reference genome
			if(referenceChromosomes.size() > 1) {
				species2errorRate.put("Reference", new HashMap<Integer,Double>());
				species2trivialName.put("Reference", "Reference");
			}
			
			
			//now generate interval trees on all contained genomes
			HashMap<String,IntervalTree<Microbe>> index2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
			File[] indexFiles = new File(genomeIndexDirPath).listFiles();
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					IntervalTree<Microbe> genomeTree = new IntervalTree<Microbe>();
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						splittedId = barPattern.split(splittedLine[0]);
						
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
						if(species2errorRate.containsKey(shortId) || microbialContentFilePath == null) {
							Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
							genomeTree.add(tmpMicrobe);
						}
						else if(shortId.contains("NZ_")) {
							String tmpMicrobeId = shortId.substring(0,9) + "000000";
							if(species2errorRate.containsKey(tmpMicrobeId) || microbialContentFilePath == null) {
								Microbe tmpMicrobe = new Microbe(tmpMicrobeId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
								genomeTree.add(tmpMicrobe);
							}
						}
						else if(shortId.length() >= 6 && (species2errorRate.containsKey(shortId.substring(0,6) + "000000") || microbialContentFilePath == null)) {
							String tmpMicrobeId = shortId.substring(0,6) + "000000";
							Microbe tmpMicrobe = new Microbe(tmpMicrobeId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
							genomeTree.add(tmpMicrobe);
						}
					}
					br.close();
					index2genomeTree.put(indexFile.getName().substring(0, indexFile.getName().lastIndexOf('.')), genomeTree);
				}
			}
			
			//adding chromosomes of not indexed species
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			String chr;
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				if(!index2genomeTree.containsKey(chr) && (referenceChromosomes.size() < 2 || !referenceChromosomes.contains(chr))) {
					species2errorRate.put(chr, new HashMap<Integer,Double>());
					species2trivialName.put(chr,chr);
				}
			}
			br.close();
			
			//now go through the sam file and collect relevant information
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			int start;
			IntervalTree<Microbe> currentGenomeTree;
			Collection<Microbe> currentMicrobes;
			Microbe currentMicrobe;
			int mismatches = Integer.MAX_VALUE;
			int maxMismatches = Integer.MIN_VALUE;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				if(index2genomeTree.containsKey(chr)) {
					currentGenomeTree = index2genomeTree.get(chr);
					start = Integer.valueOf(splittedLine[3]);
					currentMicrobes = currentGenomeTree.getIntervalsSpanning(start, new ArrayList<Microbe>());
					if(!currentMicrobes.isEmpty()) {
						currentMicrobe = currentMicrobes.iterator().next();
						
						//additional sam fields start at index 11
						mismatches = Integer.MAX_VALUE;
						for(int i = 11; i < splittedLine.length; i++) {
							if(!useMDflag && splittedLine[i].contains("NM:i:")) {
								mismatches = Integer.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								break;
							}
							else if(useMDflag && splittedLine[i].contains("MD:Z:")) {
								mismatches = getMismatchCountFromMDtag(splittedLine[i].substring(5));
								break;
							}
						}
						
						
						if(mismatches > maxMismatches)
							maxMismatches = mismatches;
						
						
						if(microbialContentFilePath == null) {
							if(!species2errorRate.containsKey(currentMicrobe.getId())) {
								species2errorRate.put(currentMicrobe.getId(), new HashMap<Integer,Double>());
								species2trivialName.put(currentMicrobe.getId(),currentMicrobe.getId());
							}
						}
						
						if(species2errorRate.get(currentMicrobe.getId()).containsKey(mismatches))
							species2errorRate.get(currentMicrobe.getId()).put(mismatches, species2errorRate.get(currentMicrobe.getId()).get(mismatches) + 1);
						else
							species2errorRate.get(currentMicrobe.getId()).put(mismatches, 1.0);
					
					}
					
				}
				
				//found mapping on reference or to something unknown....
				else {
					mismatches = Integer.MAX_VALUE;
					for(int i = 11; i < splittedLine.length; i++) {
						if(splittedLine[i].contains("NM:i:")) {
							mismatches = Integer.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							break;
						}
					}
					
					
					if(mismatches > maxMismatches)
						maxMismatches = mismatches;
					
					tmpChr = chr;
					if(referenceChromosomes.contains(chr)) {
						tmpChr = "Reference";
					}
					
					if(species2errorRate.get(tmpChr).containsKey(mismatches))
						species2errorRate.get(tmpChr).put(mismatches, species2errorRate.get(tmpChr).get(mismatches) + 1);
					else
						species2errorRate.get(tmpChr).put(mismatches, 1.0);
				
				}
			}
			br.close();
			
			//merge yeast content (if available...)
			HashMap<Integer,Double> yeastMergedRates = new HashMap<Integer,Double>();
			int yeastCounter = 0;
			for(String speciesName : species2errorRate.keySet()) {
				if(yeastChr.contains(speciesName)) {
					yeastCounter++;
					HashMap<Integer,Double> tmpMap = species2errorRate.get(speciesName);
					for(int tmpMismatches : tmpMap.keySet()) {
						if(yeastMergedRates.containsKey(tmpMismatches)) {
							yeastMergedRates.put(tmpMismatches,yeastMergedRates.get(tmpMismatches) + tmpMap.get(tmpMismatches));
						}
						else
							yeastMergedRates.put(tmpMismatches,tmpMap.get(tmpMismatches));
					}
				}
			}
			
			if(yeastCounter > 0) {
				for(String chrName : yeastChr) {
					species2errorRate.remove(chrName);
					species2trivialName.remove(chrName);
				}
				species2errorRate.put("Yeast", yeastMergedRates);
				species2trivialName.put("Yeast","Yeast");
			}
			
			if(id2trivialName == null) {
				System.out.print("species\ttrivialname");
				for(int i = 0; i <= maxMismatches; i++)
					System.out.print("\t" + i + " mm");
				System.out.print("\n");
			
			
				for(String speciesName : species2errorRate.keySet()) {
					System.out.print(speciesName + "\t" + species2trivialName.get(speciesName));
					for(int i = 0; i <= maxMismatches; i++) {
						if(species2errorRate.get(speciesName).containsKey(i))
							System.out.print("\t" + species2errorRate.get(speciesName).get(i));
						else
							System.out.print("\t0");
					}
					System.out.print("\n");
				}
			}
			return species2errorRate;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private static int getMismatchCountFromMDtag(String MDtag) {
		int mismatches = 0;
		HashSet<String> knownCharacters = new HashSet<String>();
		knownCharacters.add("A");
		knownCharacters.add("C");
		knownCharacters.add("G");
		knownCharacters.add("T");
		knownCharacters.add("N");
		knownCharacters.add("^");
		
		boolean indelMode = false;
		
		for(int i = 0; i < MDtag.length(); i++) {
			if(knownCharacters.contains(MDtag.substring(i,i+1))) {
				if(MDtag.charAt(i) == '^')
					indelMode = true;
				
				if(!indelMode)
					mismatches++;
			}
			else
				indelMode = false;
		}
		
		return mismatches;
	}
	
	
	
	private static void averageErrorRateTables(ArrayList<String> filePaths, int maxMismatches) {
		try {
			HashMap<String,Microbe> id2microbe = new HashMap<String,Microbe>();
			HashMap<String,String> id2trivialName = new HashMap<String,String>();
			BufferedReader br;
			String currentLine;
			Pattern tabPattern = Pattern.compile("\t");
			String[] splittedLine;
			String speciesId;
			String trivialName;
			double coverage;
			Microbe currentMicrobe;
			int fileCount = 0;
			ArrayList<Integer> mismatches = new ArrayList<Integer>();
			String header = null;
			for(String filePath : filePaths) {
				
				if(!new File(filePath).exists()) {
					System.err.println("Warning file did not exist: " + filePath);
					continue;
				}
				fileCount++;
				br = new BufferedReader(new FileReader(new File(filePath)));
				//skip the header
				header = br.readLine();
				while((currentLine = br.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					speciesId = splittedLine[0];
					trivialName = splittedLine[1];
					if(!id2trivialName.containsKey(speciesId))
						id2trivialName.put(speciesId,trivialName);
					
					coverage = Double.valueOf(splittedLine[2]);
					mismatches.clear();
					for(int i = 3; i <= 3 + maxMismatches; i++) {
						mismatches.add(Integer.valueOf(splittedLine[i]));
					}
					
					if(id2microbe.containsKey(speciesId)) {
						currentMicrobe = id2microbe.get(speciesId);
						currentMicrobe.incrementCoverage(coverage);
						for(int i = 0; i < mismatches.size(); i++) {
							currentMicrobe.addMismatch(i, mismatches.get(i));
						}
					}
					else {
						currentMicrobe = new Microbe(speciesId,0,0,maxMismatches);
						currentMicrobe.setCoverage(coverage);
						for(int i = 0; i < mismatches.size(); i++) {
							currentMicrobe.addMismatch(i, mismatches.get(i));
						}
						id2microbe.put(speciesId, currentMicrobe);
					}
				}
				br.close();
			}
			
			//print merged error table
			System.err.println(fileCount);
			System.out.println(header);
			int[] currentMismatches;
			for(Microbe microbe : id2microbe.values()) {
				System.out.print(String.format("%s\t%s\t%s",microbe.getId(),id2trivialName.get(microbe.getId()),microbe.getCoverage()/(double)fileCount));
				currentMismatches = microbe.getMismatches();
				for(int i = 0; i < currentMismatches.length; i++) {
					System.out.print("\t" + currentMismatches[i]);
				}
				System.out.print("\n");
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public static void generateDistanceTableFromConcatenatedGenomes(String samFilePath, String microbialContentFilePath, String genomeIndexDirPath,HashSet<String> speciesFilter,boolean mergeYeastChr) {
		try {
			ArrayList<String> speciesList = new ArrayList<String>();
			HashMap<String,String> species2trivialName = new HashMap<String,String>();
			
			BufferedReader br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
			//skipping the header
			br.readLine();
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				
				if(mergeYeastChr && yeastChr.contains(splittedLine[0])) {
					if(!speciesList.contains("Yeast")) {
						speciesList.add("Yeast");
						species2trivialName.put("Yeast","Yeast");
					}
				}
				else {
					speciesList.add(splittedLine[0]);
					species2trivialName.put(splittedLine[0], splittedLine[1]);
				}
			}
			br.close();
			
			// first field of double array contains the number of second best hits and the second field the sum of score differences
			HashMap<String,HashMap<String,Double[]>> species2secondBestSpecies = new HashMap<String,HashMap<String,Double[]>>();
			HashMap<String,Double> species2bestHitCounter = new HashMap<String,Double>();
			for(int i = 0; i < speciesList.size(); i++) {
				species2secondBestSpecies.put(speciesList.get(i), new HashMap<String,Double[]>());
				species2bestHitCounter.put(speciesList.get(i),0.0);
				for(int j = 0; j < speciesList.size(); j++) {
					if(i != j) {
						Double[] tmpArray = {0.0,0.0};
						species2secondBestSpecies.get(speciesList.get(i)).put(speciesList.get(j), tmpArray);
					}
				}
			}
			
			
			//now generate interval trees on all contained genomes
			HashMap<String,IntervalTree<Microbe>> index2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
			File[] indexFiles = new File(genomeIndexDirPath).listFiles();
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			String shortId;
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					IntervalTree<Microbe> genomeTree = new IntervalTree<Microbe>();
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						splittedId = barPattern.split(splittedLine[0]);
						
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
					
						if(shortId.contains("NZ_")) {
							String tmpMicrobeId = shortId.substring(0,9) + "000000";
							genomeTree.add(new Microbe(tmpMicrobeId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0));
						}
						else if(shortId.length() >= 6 && species2secondBestSpecies.containsKey(shortId.substring(0,6) + "000000")) {
							String tmpMicrobeId = shortId.substring(0,6) + "000000";
							genomeTree.add(new Microbe(tmpMicrobeId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0));
						}
						
						else if(mergeYeastChr && yeastChr.contains(shortId)) {
							genomeTree.add(new Microbe("Yeast",Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0));
						}
						else {
							Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
							genomeTree.add(tmpMicrobe);
						}
					}
					br.close();
					index2genomeTree.put(indexFile.getName().substring(0, indexFile.getName().lastIndexOf('.')), genomeTree);
				}
			}
						
			//now go through the sam file and collect relevant information
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			String chr;
			int start;
			IntervalTree<Microbe> currentGenomeTree;
			Collection<Microbe> currentMicrobes;
			Microbe bestHitMicrobe;
			Microbe secondBestHitMicrobe;
			String chrOfSecondBestHit;
			int startOfSecondBestHit = -1;
			double scoreOfBestHit = Double.MIN_VALUE;
			double scoreOfSecondBestHit = Double.MIN_VALUE;
			double scoreDifference;
			HashMap<String,Double[]> currentSpeciesMap;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				if(index2genomeTree.containsKey(chr)) {
					currentGenomeTree = index2genomeTree.get(chr);
					start = Integer.valueOf(splittedLine[3]);
					currentMicrobes = currentGenomeTree.getIntervalsSpanning(start, new ArrayList<Microbe>());
					
					if(!currentMicrobes.isEmpty()) {
						bestHitMicrobe = currentMicrobes.iterator().next();
						currentSpeciesMap = species2secondBestSpecies.get(bestHitMicrobe.getId());
						species2bestHitCounter.put(bestHitMicrobe.getId(),species2bestHitCounter.get(bestHitMicrobe.getId()) + 1.0);
						
						
						chrOfSecondBestHit = null;
						for(int i = 11; i < splittedLine.length; i++) {
							if(splittedLine[i].contains("CC:Z:")) {
								chrOfSecondBestHit = splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1);
								continue;
							}
							if(splittedLine[i].contains("CP:i:")) {
								startOfSecondBestHit = Integer.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								continue;
							}
							if(splittedLine[i].contains("S1:f:")) {
								scoreOfBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								continue;
							}
							if(splittedLine[i].contains("S2:f:")) {
								scoreOfSecondBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								continue;
							}
							
						}
						if(chrOfSecondBestHit != null && index2genomeTree.containsKey(chrOfSecondBestHit)) {
							currentGenomeTree = index2genomeTree.get(chrOfSecondBestHit);
							currentMicrobes = currentGenomeTree.getIntervalsSpanning(startOfSecondBestHit, new ArrayList<Microbe>());
							if(!currentMicrobes.isEmpty()) {
								secondBestHitMicrobe = currentMicrobes.iterator().next();
								if(!bestHitMicrobe.getId().equals(secondBestHitMicrobe.getId()) && currentSpeciesMap.containsKey(secondBestHitMicrobe.getId())) {
									if(scoreOfBestHit != 0.0) {
										scoreDifference = Math.abs(scoreOfBestHit - scoreOfSecondBestHit)/Math.abs(scoreOfBestHit);
									}
									else
										scoreDifference = 0.0;
									
									currentSpeciesMap.get(secondBestHitMicrobe.getId())[0]++;
									currentSpeciesMap.get(secondBestHitMicrobe.getId())[1]+=scoreDifference;
								}
								
							}
							else {
								System.err.println(String.format("WARNING: Found no entry for chr '%s' at position '%s'",chr,startOfSecondBestHit));
							}
						}
					}
					else {
						System.err.println(String.format("WARNING: Found no entry for chr '%s' at position '%s'",chr,start));
					}
					
				}
			}
			br.close();
			
			
			//output distance table now
			System.out.print("Distance");
			for(int i = 0; i < speciesList.size(); i++) {
				System.out.print("\t" + speciesList.get(i) + "|" + species2trivialName.get(speciesList.get(i)));
			}
			System.out.print("\n");
			
			
			double distance;
			double bestHitCounter;
			double secondBestHitCounter;
			double sumOfScoreDifferences;
			
			double testCounter = 0;
			for(int i = 0;i < speciesList.size(); i++) {
				if(speciesFilter.size() > 0 && !speciesFilter.contains(speciesList.get(i)))
					continue;
				
				System.out.print(speciesList.get(i) + "|" + species2trivialName.get(speciesList.get(i)));
				
				for(int j = 0; j < speciesList.size(); j++) {
					if(i == j)
						distance = 0.0;
					else {
						
						//calculate first part of distance
						
						bestHitCounter = species2bestHitCounter.get(speciesList.get(i));
						secondBestHitCounter = species2secondBestSpecies.get(speciesList.get(i)).get(speciesList.get(j))[0];
						
						sumOfScoreDifferences = species2secondBestSpecies.get(speciesList.get(i)).get(speciesList.get(j))[1];
						if(secondBestHitCounter != 0) {
							distance = (sumOfScoreDifferences + bestHitCounter - secondBestHitCounter)/bestHitCounter;
							//distance = (1 - (secondBestHitCounter/bestHitCounter)) * (sumOfScoreDifferences/(secondBestHitCounter));
						}
						else
							distance = 1.0;
						
						//calculate second part of distance
						bestHitCounter = species2bestHitCounter.get(speciesList.get(j));
						secondBestHitCounter = species2secondBestSpecies.get(speciesList.get(j)).get(speciesList.get(i))[0];
						sumOfScoreDifferences = species2secondBestSpecies.get(speciesList.get(j)).get(speciesList.get(i))[1];
						if(secondBestHitCounter != 0) {
							distance += (sumOfScoreDifferences + bestHitCounter - secondBestHitCounter)/bestHitCounter;
							//distance += (1 - (secondBestHitCounter/bestHitCounter)) * (sumOfScoreDifferences/(secondBestHitCounter));
							
						}
						else
							distance += 1.0;
							
						distance = distance/2;
						
					}
					System.out.print("\t" + trimDouble(distance));
					
				}
				System.out.print("\n");
			}
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void generateDistanceTableFromUnconcatenatedGenomes(String samFilePath, String microbialContentFilePath, String genomeIndexDirPath,HashSet<String> speciesFilter,boolean mergeYeastChr) {
		try {
			ArrayList<String> speciesList = new ArrayList<String>();
			HashMap<String,String> species2trivialName = new HashMap<String,String>();
			
			BufferedReader br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
			//skipping the header
			br.readLine();
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				
				if(mergeYeastChr && yeastChr.contains(splittedLine[0])) {
					if(!speciesList.contains("Yeast")) {
						speciesList.add("Yeast");
						species2trivialName.put("Yeast","Yeast");
					}
				}
				else {
					speciesList.add(splittedLine[0]);
					species2trivialName.put(splittedLine[0], splittedLine[1]);
				}
			}
			br.close();
			
			// first field of double array contains the number of second best hits and the second field the sum of score differences
			HashMap<String,HashMap<String,Double[]>> species2secondBestSpecies = new HashMap<String,HashMap<String,Double[]>>();
			HashMap<String,Double> species2bestHitCounter = new HashMap<String,Double>();
			for(int i = 0; i < speciesList.size(); i++) {
				species2secondBestSpecies.put(speciesList.get(i), new HashMap<String,Double[]>());
				species2bestHitCounter.put(speciesList.get(i),0.0);
				for(int j = 0; j < speciesList.size(); j++) {
					if(i != j) {
						Double[] tmpArray = {0.0,0.0};
						species2secondBestSpecies.get(speciesList.get(i)).put(speciesList.get(j), tmpArray);
					}
				}
			}
			
			
			//hash all indexed microbes
			HashMap<String,Microbe> idx2microbe = new HashMap<String,Microbe>();
			File[] indexFiles = new File(genomeIndexDirPath).listFiles();
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			String indexId;
			String shortId;
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						indexId = splittedLine[0];
						splittedId = barPattern.split(splittedLine[1]);
						
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
						Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[2]),0);
						idx2microbe.put(indexId,tmpMicrobe);
					}
					br.close();
				}
			}
						
			//now go through the sam file and collect relevant information
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			String chr;
			int start = -1;
			IntervalTree<Microbe> currentGenomeTree;
			Collection<Microbe> currentMicrobes;
			Microbe bestHitMicrobe;
			Microbe secondBestHitMicrobe;
			String chrOfSecondBestHit;
			int startOfSecondBestHit = -1;
			double scoreOfBestHit = Double.MIN_VALUE;
			double scoreOfSecondBestHit = Double.MIN_VALUE;
			double scoreDifference;
			HashMap<String,Double[]> currentSpeciesMap;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				if(idx2microbe.containsKey(chr)) {
					bestHitMicrobe = idx2microbe.get(chr);
					currentSpeciesMap = species2secondBestSpecies.get(bestHitMicrobe.getId());
					species2bestHitCounter.put(bestHitMicrobe.getId(),species2bestHitCounter.get(bestHitMicrobe.getId()) + 1.0);
					
					
					chrOfSecondBestHit = null;
					for(int i = 11; i < splittedLine.length; i++) {
						if(splittedLine[i].contains("CC:Z:")) {
							chrOfSecondBestHit = splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1);
							continue;
						}
						if(splittedLine[i].contains("CP:i:")) {
							startOfSecondBestHit = Integer.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							continue;
						}
						if(splittedLine[i].contains("S1:f:")) {
							scoreOfBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							continue;
						}
						if(splittedLine[i].contains("S2:f:")) {
							scoreOfSecondBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							continue;
						}
						
					}
					if(chrOfSecondBestHit != null && idx2microbe.containsKey(chrOfSecondBestHit)) {
						secondBestHitMicrobe = idx2microbe.get(chrOfSecondBestHit);
						if(!bestHitMicrobe.getId().equals(secondBestHitMicrobe.getId()) && currentSpeciesMap.containsKey(secondBestHitMicrobe.getId())) {
							if(scoreOfBestHit != 0.0) {
								scoreDifference = Math.abs(scoreOfBestHit - scoreOfSecondBestHit)/Math.abs(scoreOfBestHit);
							}
							else
								scoreDifference = 0.0;
							
							currentSpeciesMap.get(secondBestHitMicrobe.getId())[0]++;
							currentSpeciesMap.get(secondBestHitMicrobe.getId())[1]+=scoreDifference;
						}
					}
				
					else if(chrOfSecondBestHit != null) {
						System.err.println(String.format("WARNING: Found no entry for chr '%s' at position '%s'",chr,start));
					}
					
				
				}
			}
			br.close();
			
			
			//output distance table now
			System.out.print("Distance");
			for(int i = 0; i < speciesList.size(); i++) {
				System.out.print("\t" + speciesList.get(i) + "|" + species2trivialName.get(speciesList.get(i)));
			}
			System.out.print("\n");
			
			
			double distance;
			double bestHitCounter;
			double secondBestHitCounter;
			double sumOfScoreDifferences;
			
			double testCounter = 0;
			for(int i = 0;i < speciesList.size(); i++) {
				if(speciesFilter.size() > 0 && !speciesFilter.contains(speciesList.get(i)))
					continue;
				
				System.out.print(speciesList.get(i) + "|" + species2trivialName.get(speciesList.get(i)));
				
				for(int j = 0; j < speciesList.size(); j++) {
					if(i == j)
						distance = 0.0;
					else {
						
						//calculate first part of distance
						
						bestHitCounter = species2bestHitCounter.get(speciesList.get(i));
						secondBestHitCounter = species2secondBestSpecies.get(speciesList.get(i)).get(speciesList.get(j))[0];
						
						sumOfScoreDifferences = species2secondBestSpecies.get(speciesList.get(i)).get(speciesList.get(j))[1];
						if(secondBestHitCounter != 0) {
							distance = (sumOfScoreDifferences + bestHitCounter - secondBestHitCounter)/bestHitCounter;
							//distance = (1 - (secondBestHitCounter/bestHitCounter)) * (sumOfScoreDifferences/(secondBestHitCounter));
						}
						else
							distance = 1.0;
						
						//calculate second part of distance
						bestHitCounter = species2bestHitCounter.get(speciesList.get(j));
						secondBestHitCounter = species2secondBestSpecies.get(speciesList.get(j)).get(speciesList.get(i))[0];
						sumOfScoreDifferences = species2secondBestSpecies.get(speciesList.get(j)).get(speciesList.get(i))[1];
						if(secondBestHitCounter != 0) {
							distance += (sumOfScoreDifferences + bestHitCounter - secondBestHitCounter)/bestHitCounter;
							//distance += (1 - (secondBestHitCounter/bestHitCounter)) * (sumOfScoreDifferences/(secondBestHitCounter));
							
						}
						else
							distance += 1.0;
							
						distance = distance/2;
						
					}
					System.out.print("\t" + trimDouble(distance));
					
				}
				System.out.print("\n");
			}
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void filterDistanceTable(String tableFilePath, HashSet<String> speciesNames, double maxDistance) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(tableFilePath)));
			HashSet<Integer> columnsOfInterest = new HashSet<Integer>();
			//reading the header and getting columns of interest
			String currentLine = br.readLine();
			Pattern tabPattern = Pattern.compile("\t");
			String[] splittedLine = tabPattern.split(currentLine);
			
			for(int i = 0; i < splittedLine.length; i++) {
				if(speciesNames.contains(splittedLine[i])) {
					columnsOfInterest.add(i);
				}
			}
			
			
			ArrayList<Integer> additionalColumns = new ArrayList<Integer>();
			int currentColumn = 0;
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				if(columnsOfInterest.contains(++currentColumn)) {
					for(int i = 1; i < splittedLine.length; i++) {
						if(Double.valueOf(splittedLine[i]) < maxDistance) {
							additionalColumns.add(i);
						}
					}
				}
			}
			br.close();
			columnsOfInterest.addAll(additionalColumns);
			
			//write filtered table
			ArrayList<Integer> columns = new ArrayList<Integer>();
			columns.addAll(columnsOfInterest);
			columnsOfInterest = null;
			
			Collections.sort(columns);
			br = new BufferedReader(new FileReader(new File(tableFilePath)));
			currentLine = br.readLine();
			splittedLine = tabPattern.split(currentLine);
			System.out.print("Distance");
			for(int col : columns) {
				System.out.print("\t" + splittedLine[col]);
			}
			System.out.println();
			int currentRow = 0;
			while((currentLine = br.readLine()) != null) {
				if(columns.contains(++currentRow)) {
					splittedLine = tabPattern.split(currentLine);
					System.out.print(splittedLine[0]);
					for(int col : columns) {
						System.out.print("\t" + splittedLine[col]);
					}
					System.out.println();
				}
			}
			br.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static HashMap<String,Double> generateConfidenceScoresFromConcatenatedGenomes(String samFilePath, String microbialContentFilePath, String genomeIndexDirPath, boolean mergeYeastChr, HashMap<String,String> hashedMicrobialContent) {
		try {
			Date date;
			HashMap<String,Double> species2confidence = new HashMap<String,Double>();
			BufferedReader br;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			ArrayList<String> speciesList = new ArrayList<String>();
			HashMap<String,String> species2trivialName = new HashMap<String,String>();

			if(hashedMicrobialContent != null) {
				species2trivialName = hashedMicrobialContent;
				for(String species : species2trivialName.keySet()) {
					if(mergeYeastChr && yeastChr.contains(species)) {
						if(!speciesList.contains("Yeast")) {
							speciesList.add("Yeast");
						}
					}
					else {
						speciesList.add(species);
					}
				}
				if(mergeYeastChr) {
					species2trivialName.put("Yeast","Yeast");
					for(String chr : yeastChr) {
						species2trivialName.remove(chr);
					}
				}
			}
			else {
				date = new Date();
				System.err.println(String.format("[%s]\tReading species list.",date.toLocaleString()));
				br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
				//skipping the header
				br.readLine();
				
				while(br.ready()) {
					currentLine = br.readLine();
					splittedLine = tabPattern.split(currentLine);
					
					if(mergeYeastChr && yeastChr.contains(splittedLine[0])) {
						if(!speciesList.contains("Yeast")) {
							speciesList.add("Yeast");
							species2trivialName.put("Yeast","Yeast");
						}
					}
					else {
						speciesList.add(splittedLine[0]);
						species2trivialName.put(splittedLine[0], splittedLine[1]);
					}
				}
				br.close();
			}
			// first field of double array contains the number of second best hits and the second field the sum of score differences
			HashMap<String,HashMap<String,Double[]>> species2secondBestSpecies = new HashMap<String,HashMap<String,Double[]>>();
			HashMap<String,Double> species2bestHitCounter = new HashMap<String,Double>();
			for(int i = 0; i < speciesList.size(); i++) {
				species2secondBestSpecies.put(speciesList.get(i), new HashMap<String,Double[]>());
				species2bestHitCounter.put(speciesList.get(i),0.0);
				for(int j = 0; j < speciesList.size(); j++) {
					if(i != j) {
						Double[] tmpArray = {0.0,0.0};
						species2secondBestSpecies.get(speciesList.get(i)).put(speciesList.get(j), tmpArray);
					}
				}
			}
			
			date = new Date();
			if(hashedMicrobialContent == null)
				System.err.println(String.format("[%s]\tReading index files.",date.toLocaleString()));
			//now generate interval trees on all contained genomes
			HashMap<String,IntervalTree<Microbe>> index2genomeTree = new HashMap<String,IntervalTree<Microbe>>();
			File[] indexFiles = new File(genomeIndexDirPath).listFiles();
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			String shortId;
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					IntervalTree<Microbe> genomeTree = new IntervalTree<Microbe>();
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						splittedId = barPattern.split(splittedLine[0]);
						
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
						if(shortId.contains("NZ_")) {
							String tmpMicrobeId = shortId.substring(0,9) + "000000";
							genomeTree.add(new Microbe(tmpMicrobeId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0));
						}
						else if(shortId.length() >= 6 && species2secondBestSpecies.containsKey(shortId.substring(0,6) + "000000")) {
							String tmpMicrobeId = shortId.substring(0,6) + "000000";
							genomeTree.add(new Microbe(tmpMicrobeId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0));
						}
						else if(mergeYeastChr && yeastChr.contains(shortId)) {
							genomeTree.add(new Microbe("Yeast",Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0));
						}
						else {
							Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[1]),Integer.valueOf(splittedLine[2]),0);
							genomeTree.add(tmpMicrobe);
						}
						
					}
					br.close();
					index2genomeTree.put(indexFile.getName().substring(0, indexFile.getName().lastIndexOf('.')), genomeTree);
				}
			}
			
			
			date = new Date();
			if(hashedMicrobialContent == null)
				System.err.println(String.format("[%s]\tProcessing sam file.",date.toLocaleString()));
			//now go through the sam file and collect relevant information
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			String chr;
			int start;
			IntervalTree<Microbe> currentGenomeTree;
			Collection<Microbe> currentMicrobes;
			Microbe bestHitMicrobe;
			Microbe secondBestHitMicrobe;
			String chrOfSecondBestHit;
			int startOfSecondBestHit = -1;
			double scoreOfBestHit = Double.MIN_VALUE;
			double scoreOfSecondBestHit = Double.MIN_VALUE;
			double scoreDifference;
			HashMap<String,Double[]> currentSpeciesMap;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				if(index2genomeTree.containsKey(chr)) {
					currentGenomeTree = index2genomeTree.get(chr);
					start = Integer.valueOf(splittedLine[3]);
					currentMicrobes = currentGenomeTree.getIntervalsSpanning(start, new ArrayList<Microbe>());
					if(!currentMicrobes.isEmpty()) {
						bestHitMicrobe = currentMicrobes.iterator().next();						
						currentSpeciesMap = species2secondBestSpecies.get(bestHitMicrobe.getId());
						species2bestHitCounter.put(bestHitMicrobe.getId(),species2bestHitCounter.get(bestHitMicrobe.getId()) + 1.0);
						
						
						chrOfSecondBestHit = null;
						for(int i = 11; i < splittedLine.length; i++) {
							if(splittedLine[i].contains("CC:Z:")) {
								chrOfSecondBestHit = splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1);
								continue;
							}
							if(splittedLine[i].contains("CP:i:")) {
								startOfSecondBestHit = Integer.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								continue;
							}
							if(splittedLine[i].contains("S1:f:")) {
								scoreOfBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								if(scoreOfBestHit < 0)
									scoreOfBestHit = 0;
								continue;
							}
							if(splittedLine[i].contains("S2:f:")) {
								scoreOfSecondBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
								if(scoreOfSecondBestHit < 0)
									scoreOfSecondBestHit = 0;
								continue;
							}
							
						}
						if(chrOfSecondBestHit != null && index2genomeTree.containsKey(chrOfSecondBestHit)) {
							currentGenomeTree = index2genomeTree.get(chrOfSecondBestHit);
							currentMicrobes = currentGenomeTree.getIntervalsSpanning(startOfSecondBestHit, new ArrayList<Microbe>());
							if(!currentMicrobes.isEmpty()) {
								secondBestHitMicrobe = currentMicrobes.iterator().next();
								if(!bestHitMicrobe.getId().equals(secondBestHitMicrobe.getId())) {
									if(scoreOfBestHit != 0.0) {
										scoreDifference = Math.abs(scoreOfBestHit - scoreOfSecondBestHit)/Math.abs(scoreOfBestHit);
									}
									else
										scoreDifference = 0.0;
									
									if(!currentSpeciesMap.containsKey(secondBestHitMicrobe.getId())) {
										Double[] tmpArray = {0.0,0.0};
										currentSpeciesMap.put(secondBestHitMicrobe.getId(), tmpArray);
									}
									
									currentSpeciesMap.get(secondBestHitMicrobe.getId())[0]++;
									currentSpeciesMap.get(secondBestHitMicrobe.getId())[1]+=scoreDifference;
								}
							}
						}
					}
					
				}
			}
			br.close();
			
			date = new Date();
			if(hashedMicrobialContent == null)
				System.err.println(String.format("[%s]\tPrinting confidence list.",date.toLocaleString()));
			//output confidence list now.
			double confidence;
			double bestHitCounter;
			double secondBestHitCounter;
			double sumOfScoreDifferences;
			String trivialName;
			HashMap<String,Double[]> secondBestSpecies2Stats;
			if(hashedMicrobialContent == null)
				System.out.println("species\ttrivialname\tbest_hits\tsecond_best_hits\tconfidence");
			
			for(String speciesName : species2secondBestSpecies.keySet()) {
				
				trivialName = species2trivialName.get(speciesName);
				secondBestSpecies2Stats = species2secondBestSpecies.get(speciesName);
				secondBestHitCounter = 0;
				sumOfScoreDifferences = 0;
				bestHitCounter = species2bestHitCounter.get(speciesName);
				for(Double[] stats : secondBestSpecies2Stats.values()) {
					secondBestHitCounter += stats[0];
					sumOfScoreDifferences += stats[1];
				}
				
				if(secondBestHitCounter == 0)
					confidence = 1.0;
				else
					confidence = (sumOfScoreDifferences + (bestHitCounter - secondBestHitCounter))/bestHitCounter;
					//confidence = (1 - (secondBestHitCounter/bestHitCounter)) * (sumOfScoreDifferences/secondBestHitCounter);
				
				
				if(hashedMicrobialContent == null)
					System.out.println(String.format("%s\t%s\t%s\t%s\t%s",speciesName,trivialName,bestHitCounter,secondBestHitCounter,confidence));
				
				species2confidence.put(speciesName, confidence);
			}
			date = new Date();
			if(hashedMicrobialContent == null)
				System.err.println(String.format("[%s]\tDone.",date.toLocaleString()));
			return species2confidence;
		}
		
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	
	private static void generateConfidenceScoresFromUnconcatenatedGenomes(String samFilePath, String microbialContentFilePath, String genomeIndexDirPath, boolean mergeYeastChr) {
		try {
			Date date = new Date();
			System.err.println(String.format("[%s]\tReading contamination file.",date.toLocaleString()));
			ArrayList<String> speciesList = new ArrayList<String>();
			HashMap<String,String> species2trivialName = new HashMap<String,String>();
			
			BufferedReader br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
			//skipping the header
			br.readLine();
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				
				if(mergeYeastChr && yeastChr.contains(splittedLine[0])) {
					if(!speciesList.contains("Yeast")) {
						speciesList.add("Yeast");
						species2trivialName.put("Yeast","Yeast");
					}
				}
				else {
					speciesList.add(splittedLine[0]);
					species2trivialName.put(splittedLine[0], splittedLine[1]);
				}
			}
			br.close();
			
			
			// first field of double array contains the number of second best hits and the second field the sum of score differences
			HashMap<String,HashMap<String,Double[]>> species2secondBestSpecies = new HashMap<String,HashMap<String,Double[]>>();
			HashMap<String,Double> species2bestHitCounter = new HashMap<String,Double>();
			for(int i = 0; i < speciesList.size(); i++) {
				species2secondBestSpecies.put(speciesList.get(i), new HashMap<String,Double[]>());
				species2bestHitCounter.put(speciesList.get(i),0.0);
				for(int j = 0; j < speciesList.size(); j++) {
					if(i != j) {
						Double[] tmpArray = {0.0,0.0};
						species2secondBestSpecies.get(speciesList.get(i)).put(speciesList.get(j), tmpArray);
					}
				}
			}
			
			date = new Date();
			System.err.println(String.format("[%s]\tReading index files.",date.toLocaleString()));
			//now generate interval trees on all contained genomes
			HashMap<String,Microbe> idx2microbe = new HashMap<String,Microbe>();
			File[] indexFiles = new File(genomeIndexDirPath).listFiles();
			Pattern barPattern = Pattern.compile("\\|");
			String[] splittedId;
			String indexId;
			String shortId;
			for(File indexFile : indexFiles) {
				if(indexFile.getName().contains(".idx")) {
					br = new BufferedReader(new FileReader(indexFile));
					while(br.ready()) {
						currentLine = br.readLine();
						splittedLine = tabPattern.split(currentLine);
						indexId = splittedLine[0];
						splittedId = barPattern.split(splittedLine[1]);
						
						//here we have a species, which does not come from refseq
						if(splittedId.length < 3)
							shortId = splittedId[0];
						else
							shortId = splittedId[3];
							
						Microbe tmpMicrobe = new Microbe(shortId,Integer.valueOf(splittedLine[2]),0);
						idx2microbe.put(indexId,tmpMicrobe);
					}
					br.close();
					
				}
			}
			
			
			date = new Date();
			System.err.println(String.format("[%s]\tProcessing sam file.",date.toLocaleString()));
			//now go through the sam file and collect relevant information
			br = new BufferedReader(new FileReader(new File(samFilePath)));
			String chr;
			int start;
			IntervalTree<Microbe> currentGenomeTree;
			Collection<Microbe> currentMicrobes;
			Microbe bestHitMicrobe;
			Microbe secondBestHitMicrobe;
			String chrOfSecondBestHit;
			int startOfSecondBestHit = -1;
			double scoreOfBestHit = Double.MIN_VALUE;
			double scoreOfSecondBestHit = Double.MIN_VALUE;
			double scoreDifference;
			HashMap<String,Double[]> currentSpeciesMap;
			while(br.ready()) {
				currentLine = br.readLine();
				splittedLine = tabPattern.split(currentLine);
				chr = splittedLine[2];
				if(idx2microbe.containsKey(chr)) {
					bestHitMicrobe = idx2microbe.get(chr);						
					currentSpeciesMap = species2secondBestSpecies.get(bestHitMicrobe.getId());
					species2bestHitCounter.put(bestHitMicrobe.getId(),species2bestHitCounter.get(bestHitMicrobe.getId()) + 1.0);
					
					chrOfSecondBestHit = null;
					for(int i = 11; i < splittedLine.length; i++) {
						if(splittedLine[i].contains("CC:Z:")) {
							chrOfSecondBestHit = splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1);
							continue;
						}
						if(splittedLine[i].contains("CP:i:")) {
							startOfSecondBestHit = Integer.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							continue;
						}
						if(splittedLine[i].contains("S1:f:")) {
							scoreOfBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							if(scoreOfBestHit < 0)
								scoreOfBestHit = 0;
							continue;
						}
						if(splittedLine[i].contains("S2:f:")) {
							scoreOfSecondBestHit = Double.valueOf(splittedLine[i].substring(splittedLine[i].lastIndexOf(":") + 1));
							if(scoreOfSecondBestHit < 0)
								scoreOfSecondBestHit = 0;
							continue;
						}
						
					}
					if(chrOfSecondBestHit != null && idx2microbe.containsKey(chrOfSecondBestHit)) {
						
						secondBestHitMicrobe = idx2microbe.get(chrOfSecondBestHit);
						if(!bestHitMicrobe.getId().equals(secondBestHitMicrobe.getId())) {
							if(scoreOfBestHit != 0.0) {
								scoreDifference = Math.abs(scoreOfBestHit - scoreOfSecondBestHit)/Math.abs(scoreOfBestHit);
							}
							else
								scoreDifference = 0.0;
							
							if(!currentSpeciesMap.containsKey(secondBestHitMicrobe.getId())) {
								Double[] tmpArray = {0.0,0.0};
								currentSpeciesMap.put(secondBestHitMicrobe.getId(), tmpArray);
							}
							
							currentSpeciesMap.get(secondBestHitMicrobe.getId())[0]++;
							currentSpeciesMap.get(secondBestHitMicrobe.getId())[1]+=scoreDifference;
						}
						
					}
				}
			}
			br.close();
			
			date = new Date();
			System.err.println(String.format("[%s]\tPrinting confidence list.",date.toLocaleString()));
			//output confidence list now.
			double confidence;
			double bestHitCounter;
			double secondBestHitCounter;
			double sumOfScoreDifferences;
			String trivialName;
			HashMap<String,Double[]> secondBestSpecies2Stats;
			System.out.println("species\ttrivialname\tbest_hits\tsecond_best_hits\tconfidence");
			for(String speciesName : species2secondBestSpecies.keySet()) {
				trivialName = species2trivialName.get(speciesName);
				secondBestSpecies2Stats = species2secondBestSpecies.get(speciesName);
				secondBestHitCounter = 0;
				sumOfScoreDifferences = 0;
				bestHitCounter = species2bestHitCounter.get(speciesName);
				for(Double[] stats : secondBestSpecies2Stats.values()) {
					secondBestHitCounter += stats[0];
					sumOfScoreDifferences += stats[1];
				}
				
				if(secondBestHitCounter == 0)
					confidence = 1.0;
				else
					confidence = (sumOfScoreDifferences + (bestHitCounter - secondBestHitCounter))/bestHitCounter;
					//confidence = (1 - (secondBestHitCounter/bestHitCounter)) * (sumOfScoreDifferences/secondBestHitCounter);
				
				if(confidence > 1.0)
					System.err.println("WARNING: confidence score > 1.0!");
				
				System.out.println(String.format("%s\t%s\t%s\t%s\t%s",speciesName,trivialName,bestHitCounter,secondBestHitCounter,confidence));
			}
			date = new Date();
			System.err.println(String.format("[%s]\tDone.",date.toLocaleString()));
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private static void averageConfidenceScores(ArrayList<String> filePaths) {
		try {
			HashMap<String,Microbe> id2microbe = new HashMap<String,Microbe>();
			HashMap<String,String> id2trivialName = new HashMap<String,String>();
			BufferedReader br;
			String currentLine;
			Pattern tabPattern = Pattern.compile("\t");
			String[] splittedLine;
			String speciesId;
			String trivialName;
			int bestHits;
			int secondBestHits;
			double confidence;
			Microbe currentMicrobe;
			int fileCount = 0;
			String header = null;
			for(String filePath : filePaths) {
				
				if(!new File(filePath).exists()) {
					System.err.println("Warning file did not exist: " + filePath);
					continue;
				}
				fileCount++;
				br = new BufferedReader(new FileReader(new File(filePath)));
				//skip the header
				header = br.readLine();
				while((currentLine = br.readLine()) != null) {
					splittedLine = tabPattern.split(currentLine);
					speciesId = splittedLine[0];
					trivialName = splittedLine[1];
					if(!id2trivialName.containsKey(speciesId))
						id2trivialName.put(speciesId,trivialName);
					
					bestHits = Double.valueOf(splittedLine[2]).intValue();
					secondBestHits = Double.valueOf(splittedLine[3]).intValue();
					confidence = Double.valueOf(splittedLine[4]);
					
					if(id2microbe.containsKey(speciesId)) {
						currentMicrobe = id2microbe.get(speciesId);
						currentMicrobe.incrementBestHitReads(bestHits);
						currentMicrobe.incrementSecondBestHitReads(secondBestHits);
						currentMicrobe.incrementConfidence(confidence);
					}
					else {
						currentMicrobe = new Microbe(speciesId,0,0,0);
						currentMicrobe.setBestHitReads(bestHits);
						currentMicrobe.setSecondBestHitReads(secondBestHits);
						currentMicrobe.setConfidence(confidence);
						id2microbe.put(speciesId, currentMicrobe);
					}
				}
				br.close();
			}
			
			//print merged confidence table
			System.err.println(fileCount);
			System.out.println(header);
			for(Microbe microbe : id2microbe.values()) {
				System.out.println(String.format("%s\t%s\t%s\t%s\t%s",microbe.getId(),id2trivialName.get(microbe.getId()),(double)microbe.getBestHitReads()/(double)fileCount,(double)microbe.getSecondBestHitReads()/(double)fileCount,microbe.getConfidence()/(double)fileCount));
			}
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * extracts information from the following format
	 * read 49696 gene infos
	 * got 326906 correct, multi correct 0 multi wrong 0 multi: 0
	 * processed 436183 wrong start 15434 wrong chr 4379
	 * full predicted 340472 correct 326906 onesideok 248 bothsideok 0
	 * nsplit predicted 95711 correct 81323 oneside 5754 bothside(splic error) 1095
	 * FULL: FULL: 340815 predicted: 340472 (99.90%) correct: 326906 (95.92%) onesideok: 248(0.07%) bothsideok: 0(0.00%) total ok: 327154(95.99%)
	 * SPLIT: FULL: 105583 predicted: 95711 (90.65%) correct: 81323 (77.02%) onesideok: 5754(5.45%) bothsideok: 1095(1.04%) total ok: 88172(83.51%)
	 * TOTAL: FULL: 446398 predicted: 436183 (97.71%) correct: 408229 (91.45%) onesideok: 6002(1.34%) bothsideok: 1095(0.25%) total ok: 415326(93.04%)
	 * SINGLE PPV: 415326/416370=99.75% (relaxed)
	 * 
	 * extracted fields are (percentages only):
	 * FULL: correct + total ok
	 * SPLIT: correct + total ok
	 * 
	 * Numbers are printed to std out in the following order:
	 * full correct	full total ok	split correct	split total ok
	 */
	
	private static void extractAccuracyFromTable(String tableFilePath) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(tableFilePath)));
			String currentLine;
			String splittedLine[];
			while((currentLine = br.readLine()) != null) {
				if(currentLine.substring(0,5).equals("FULL:")) {
					splittedLine = currentLine.split(" ");
					System.out.print(splittedLine[8].substring(1,splittedLine[6].length() - 2) + "\t" +
							splittedLine[15].substring(splittedLine[15].indexOf('(') + 1,splittedLine[15].lastIndexOf(')') - 1) + "\t");
				}
				if(currentLine.substring(0,6).equals("SPLIT:")) {
					splittedLine = currentLine.split(" ");
					System.out.print(splittedLine[8].substring(1,splittedLine[6].length() - 2) + "\t" +
							splittedLine[15].substring(splittedLine[15].indexOf('(') + 1,splittedLine[15].lastIndexOf(')') - 1));
				}
			}
			br.close();
		}
		catch(Exception e) {
			
		}
	}
	
	
	private static void buildWorkloadStatistic(String mpstatFilePath, String outputFilePath, int maxcores) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(mpstatFilePath)));
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			pw.print("time_point");
			for(int i = 1; i <= maxcores;i++) {
				pw.print("\tuser_load_core_" + i);
			}
			pw.print("\n");
			String currentLine;
			String[] splittedLine;
			ArrayList<Double> userTimes = new ArrayList<Double>();
			Double userTime;
			int timePoint = 0;
			while((currentLine = br.readLine()) != null) {
				if(currentLine.contains("Average:")) {
					userTimes.clear();
					//skip header and summary row
					br.readLine();
					currentLine = br.readLine();
					while(currentLine != null && currentLine.contains("Average:")) {
						userTime = Double.valueOf(currentLine.substring(currentLine.indexOf('.') -3,currentLine.indexOf('.') + 3));
						userTimes.add(userTime);
						currentLine = br.readLine();
					}
					
					Collections.sort(userTimes);
					pw.print(++timePoint);
					for(int i = 1; i <= maxcores;i++) {
						pw.print("\t" + userTimes.get(userTimes.size() - i));
					}
					pw.print("\n");
				}
			}
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	private static String trimDouble(double inValue){
		DecimalFormat twoDec = new DecimalFormat("0.00",new DecimalFormatSymbols(Locale.US));
		twoDec.setGroupingUsed(false);
		return twoDec.format(inValue);
		}
	
	private static class Microbe implements Interval {
		
		private String id;
		private int start;
		private int stop;
		private int genomeLength;
		private int bestHitReads;
		private int secondBestHitReads;
		private double confidence;
		private double coverage;
		
		private int[] mismatches;
		
		public Microbe(String id, int start, int end, int maxMismatches) {
			this.id = id;
			this.start = start;
			this.stop = end;
			this.genomeLength = stop - start + 1;
			this.bestHitReads = 0;
			this.secondBestHitReads = 0;
			this.confidence = 0;
			this.coverage = 0;
			
			this.mismatches = new int[maxMismatches + 1];
		}
		
		public Microbe(String id, int length, int maxMismatches) {
			this.id = id;
			this.genomeLength = length;
			this.bestHitReads = 0;
			this.secondBestHitReads = 0;
			this.confidence = 0;
			this.coverage = 0;
			
			this.mismatches = new int[maxMismatches + 1];
		}
		
		public String getId() {
			return this.id;
		}
		
		public int getStart() {
			return this.start;
		}
		
		public int getStop() {
			return this.stop;
		}
		
		public int getLength() {
			return this.genomeLength;
		}
		
		public int getBestHitReads() {
			return this.bestHitReads;
		}
		
		public int getSecondBestHitReads() {
			return this.secondBestHitReads;
		}
		
		public void setBestHitReads(int reads) {
			this.bestHitReads = reads;
		}
		
		public void incrementBestHitReads() {
			this.bestHitReads++;
		}
		
		public void incrementBestHitReads(int value) {
			this.bestHitReads+=value;
		}
		
		public void setSecondBestHitReads(int reads) {
			this.secondBestHitReads = reads;
		}
		
		public void incrementSecondBestHitReads() {
			this.secondBestHitReads++;
		}
		public void incrementSecondBestHitReads(int value) {
			this.secondBestHitReads+=value;
		}
		
		public double getConfidence() {
			return this.confidence;
		}
		
		
		public void setConfidence(double value) {
			this.confidence = value;
		}
		
		public void incrementConfidence(double value) {
			this.confidence += value;
		}
		
		public void setCoverage(double value) {
			this.coverage = value;
		}
		
		public void incrementCoverage(double value) {
			this.coverage += value;
		}
		
		public double getCoverage() {
			return this.coverage;
		}
		
		
		public void addMismatch(int mismatchNumber, int mismatchCount) {
			this.mismatches[mismatchNumber] += mismatchCount;
		}
		
		public int[] getMismatches() {
			return this.mismatches;
		}
		
	}
	
	
}
