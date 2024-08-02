package tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.regex.Pattern;

public class GenomeIndexer {

	private String genomesFilePath;
	private String outputDirPath;
	private int gapSize;
	private int maxBasePairsPerIndex;
	
	public GenomeIndexer(String genomesFilePath, String outputDirPath, int gapSize, int maxBasePairsPerIndex) {
		this.genomesFilePath = genomesFilePath;
		this.outputDirPath = outputDirPath;
		
		this.gapSize = gapSize;
		this.maxBasePairsPerIndex = maxBasePairsPerIndex;
	}
	
	public void indexWithoutConcatenation(String genomePrefix, long maxBasePairsPerIndex, String microbialContentFilePath) {
		try {
			HashSet<String> positiveList = new HashSet<String>();
			if(microbialContentFilePath != null) {
				BufferedReader br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
				while(br.ready()) {
					positiveList.add(br.readLine().split("\t")[1]);
				}
				br.close();
			}
			
			BufferedReader br = new BufferedReader(new FileReader(new File(this.genomesFilePath)));
			int sequenceCounter = 0;
			int fileCounter = 0;
			PrintWriter currentSequencePw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/%s_%s.fa",genomePrefix,fileCounter))));
			PrintWriter currentIndexPw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/%s.idx",genomePrefix))));
			long baseCounter = 0;
			long genomeLengthCounter = 0;
			String currentLine = br.readLine();
			Pattern spacePattern = Pattern.compile(" ");
			String[] splittedLine = spacePattern.split(currentLine);
			String currentSpecies = splittedLine[0].substring(1);
			StringBuilder trivialNameBuilder = new StringBuilder();
			trivialNameBuilder.append(splittedLine[1]);
			String currentTrivialName;
			for(int i = 2; i < splittedLine.length;i++) {
				if(splittedLine[i].charAt(splittedLine[i].length()-1) == ',') {
					trivialNameBuilder.append(" ").append(splittedLine[i].substring(0,splittedLine[i].length()-1));
					break;
				}
				trivialNameBuilder.append(" ").append(splittedLine[i]);
			}
			currentTrivialName = trivialNameBuilder.toString();
			boolean writeSequence = true;
			if(microbialContentFilePath != null && !positiveList.contains(currentTrivialName))
				writeSequence = false;
			
			if(writeSequence)
				currentSequencePw.println(String.format(">%s_%s",genomePrefix,sequenceCounter));
			while(br.ready()) {
				currentLine = br.readLine();
				//write current species here
				if(currentLine.charAt(0) == '>') {
					//first write the index of the previous header
					if(writeSequence)
						currentIndexPw.println(String.format("%s_%s\t%s\t%s\t%s",genomePrefix,sequenceCounter,currentSpecies,genomeLengthCounter,currentTrivialName));
					
					//start a new index
					if(baseCounter > maxBasePairsPerIndex) {
						currentSequencePw.close();
						baseCounter = 0;
						fileCounter++;
						currentSequencePw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/%s_%s.fa",genomePrefix,fileCounter))));
					}
					
					
					splittedLine = spacePattern.split(currentLine);
					currentSpecies = splittedLine[0].substring(1);
					trivialNameBuilder.setLength(0);
					trivialNameBuilder.append(splittedLine[1]);
					for(int i = 2; i < splittedLine.length;i++) {
						if(splittedLine[i].charAt(splittedLine[i].length()-1) == ',') {
							trivialNameBuilder.append(" ").append(splittedLine[i].substring(0,splittedLine[i].length()-1));
							break;
						}
						trivialNameBuilder.append(" ").append(splittedLine[i]);
					}
					currentTrivialName = trivialNameBuilder.toString();
					writeSequence = true;
					if(microbialContentFilePath != null && !positiveList.contains(currentTrivialName))
						writeSequence = false;
					
					sequenceCounter++;
					if(writeSequence)
						currentSequencePw.println(String.format(">%s_%s",genomePrefix,sequenceCounter));
					genomeLengthCounter = 0;
				}
				else {
					baseCounter += currentLine.length();
					genomeLengthCounter += currentLine.length();
					if(writeSequence)
						currentSequencePw.println(currentLine);
				}
			}
			br.close();
			currentSequencePw.close();
			
			//write last index
			if(writeSequence)
				currentIndexPw.println(String.format("%s_%s\t%s\t%s\t%s",genomePrefix,sequenceCounter,currentSpecies,genomeLengthCounter,currentTrivialName));
			currentIndexPw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void indexWithConcatenation(String genomePrefix, String microbialContentFilePath) {
		try {
			HashSet<String> positiveList = new HashSet<String>();
			if(microbialContentFilePath != null) {
				BufferedReader br = new BufferedReader(new FileReader(new File(microbialContentFilePath)));
				while(br.ready()) {
					positiveList.add(br.readLine().split("\t")[1]);
				}
				br.close();
			}
			
			BufferedReader br = new BufferedReader(new FileReader(new File(this.genomesFilePath)));
			int indexCounter = 0;
			
			File tmpFile = new File(this.outputDirPath + "/fasta");
			if(!tmpFile.isDirectory())
				tmpFile.mkdirs();
			
			tmpFile = new File(this.outputDirPath + "/indices");
			if(!tmpFile.isDirectory())
				tmpFile.mkdirs();
			
			PrintWriter currentSequencePw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/fasta/%s_%s.fa",genomePrefix,indexCounter))));
			PrintWriter currentIndexPw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/indices/%s_%s.idx",genomePrefix,indexCounter))));
			int baseCounter = 0;
			int currentStart = baseCounter;
			String currentLine = br.readLine();
			Pattern spacePattern = Pattern.compile(" ");
			String[] splittedLine = spacePattern.split(currentLine);
			String currentSpecies = splittedLine[0].substring(1);
			StringBuilder trivialNameBuilder = new StringBuilder();
			trivialNameBuilder.append(splittedLine[1]);
			String currentTrivialName;
			for(int i = 2; i < splittedLine.length;i++) {
				if(splittedLine[i].charAt(splittedLine[i].length()-1) == ',') {
					trivialNameBuilder.append(" ").append(splittedLine[i].substring(0,splittedLine[i].length()-1));
					break;
				}
				trivialNameBuilder.append(" ").append(splittedLine[i]);
			}
			currentTrivialName = trivialNameBuilder.toString();
			boolean writeSequence = true;
			if(microbialContentFilePath != null && !positiveList.contains(currentTrivialName))
				writeSequence = false;
			
			StringBuilder gap = new StringBuilder();
			for(int i = 0; i < this.gapSize; i++) {
				gap.append("N");
			}
			
			currentSequencePw.println(String.format(">%s_%s",genomePrefix,indexCounter));
			while(br.ready()) {
				currentLine = br.readLine();
				//write current species here
				if(currentLine.charAt(0) == '>') {
					if(writeSequence) {
						currentIndexPw.println(String.format("%s\t%s\t%s\t%s",currentSpecies,currentStart,currentStart + baseCounter,currentTrivialName));
						//add 'N's to seperate the current species from the next
						currentSequencePw.println(gap.toString());
						baseCounter += gap.length();
					}
					//start a new index
					if(currentStart + baseCounter > this.maxBasePairsPerIndex) {
						currentSequencePw.close();
						currentIndexPw.close();
						currentStart = 0;
						
						indexCounter++;
						currentSequencePw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/fasta/%s_%s.fa",genomePrefix,indexCounter))));
						currentIndexPw = new PrintWriter(new FileWriter(new File(this.outputDirPath + String.format("/indices/%s_%s.idx",genomePrefix,indexCounter))));
						currentSequencePw.println(String.format(">%s_%s",genomePrefix,indexCounter));
					}
					
					else if(writeSequence)
						currentStart += baseCounter;
					
					splittedLine = spacePattern.split(currentLine);
					currentSpecies = splittedLine[0].substring(1);
					trivialNameBuilder.setLength(0);
					trivialNameBuilder.append(splittedLine[1]);
					for(int i = 2; i < splittedLine.length;i++) {
						if(splittedLine[i].charAt(splittedLine[i].length()-1) == ',') {
							trivialNameBuilder.append(" ").append(splittedLine[i].substring(0,splittedLine[i].length()-1));
							break;
						}
						trivialNameBuilder.append(" ").append(splittedLine[i]);
					}
					currentTrivialName = trivialNameBuilder.toString();
					writeSequence = true;
					if(microbialContentFilePath != null && !positiveList.contains(currentTrivialName))
						writeSequence = false;
					
					baseCounter = 0;
				}
				else {
					if(writeSequence) {
						baseCounter += currentLine.length();
						currentSequencePw.println(currentLine);
					}
				}
			}
			br.close();
			
			//write the last species to the index (sequence already written!)
			if(writeSequence)
				currentIndexPw.println(String.format("%s\t%s\t%s\t%s",currentSpecies,currentStart,(currentStart + baseCounter),currentTrivialName));
			
			
			currentSequencePw.close();
			currentIndexPw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
}
