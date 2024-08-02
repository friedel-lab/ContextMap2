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
import java.util.TreeMap;

import main.MutableInt;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.SeekableBufferedStream;
import net.sf.samtools.util.SeekableFileStream;

public class BamConverter {

	
	public BamConverter() {
		
	}
	
	
	public void getSeedMismatchStats(String bamFilePath, String indexFilePath, int seedLength, int maxMismatches, String outputFilePath) {
		try {
			SAMFileReader bamReader = new SAMFileReader(new SeekableBufferedStream(new SeekableFileStream(new File(bamFilePath))),new File(indexFilePath),true);
			SAMRecordIterator recordIterator = bamReader.iterator();
			SAMRecord samRecord;
			Cigar cigar;
			CigarElement currentElement;
			CigarOperator currentOperator;
			int currentLength;
			int currentMismatchCount;
			TreeMap<Integer,MutableInt> mismatch2count = new TreeMap<Integer,MutableInt>();
			PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFilePath)));
			ArrayList<CigarElement> currentElements = new ArrayList<CigarElement>();
			while(recordIterator.hasNext()) {
				samRecord = recordIterator.next();
				cigar = samRecord.getCigar();
				currentLength = 0;
				currentMismatchCount = 0;
				
				currentElements.clear();
				currentElements.addAll(cigar.getCigarElements());
				if(samRecord.getReadNegativeStrandFlag())
					Collections.reverse(currentElements);
				
				for(int i = 0; i < currentElements.size(); i++) {
					currentElement = currentElements.get(i);
					currentOperator = currentElement.getOperator();
					
					
					if(currentOperator.equals(CigarOperator.EQ)) {
						currentLength += currentElement.getLength();
					}
					
					//check if we already exceeded the seed region
					if(currentLength > seedLength) {
						break;
					}
					
					if(currentOperator.equals(CigarOperator.X)) {
						currentMismatchCount += currentElement.getLength();
					}
				}
				
				if(currentMismatchCount > maxMismatches) {
					if(mismatch2count.containsKey(currentMismatchCount)) {
						mismatch2count.get(currentMismatchCount).increment();
					}
					else {
						mismatch2count.put(currentMismatchCount,new MutableInt(1));
					}
					
					
					
					pw.println(samRecord.getReadName());
				}
				
			}
			recordIterator.close();
			pw.close();
			System.out.println("mismatches\tcount");
			for(int mismatchNumber : mismatch2count.keySet()) {
				System.out.println(String.format("%s\t%s",mismatchNumber,mismatch2count.get(mismatchNumber)));
			}
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void addMismatchPositionsToCigarString(String bamFilePath, String indexFilePath, String genomeDirPath, String fastaFilePath, String outputFilePath,boolean checkModification,boolean prebuffer) {
		try {
			SAMFileReader bamReader = new SAMFileReader(new SeekableBufferedStream(new SeekableFileStream(new File(bamFilePath))),new File(indexFilePath),true);
			SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(bamReader.getFileHeader(),true,new File(outputFilePath));
			File[] chrFiles = new File(genomeDirPath).listFiles();
			String currentChrName;
			StringBuffer chrSequence = new StringBuffer();
			StringBuffer readSequence = new StringBuffer();
			HashMap<String,String> readid2sequence = new HashMap<String,String>();
			BufferedReader sequenceReader;
			String currentLine;
			
			//prebuffering read sequences
			Date date;
			if(prebuffer) {
				date = new Date();
				System.out.println(String.format("[%s]\tPrebuffering read sequences.",date.toLocaleString()));
				sequenceReader = new BufferedReader(new FileReader(new File(fastaFilePath)));
				while((currentLine = sequenceReader.readLine()) != null) {
					if(currentLine.charAt(0) == '>')
						readid2sequence.put(currentLine.substring(1), sequenceReader.readLine());
				}
				sequenceReader.close();
			}
			
			
			SAMRecordIterator recordIterator;
			
			for(File chrFile : chrFiles) {
				//get chromosome sequence
				date = new Date();
				currentChrName = chrFile.getName().substring(0,chrFile.getName().lastIndexOf('.'));
				System.out.println(String.format("[%s]\tBuffering next chromosome sequence (%s).",date.toLocaleString(),currentChrName));
				chrSequence.setLength(0);
				sequenceReader = new BufferedReader(new FileReader(chrFile));
				sequenceReader.readLine();
				while(sequenceReader.ready()) {
					chrSequence.append(sequenceReader.readLine().toUpperCase());
				}
				sequenceReader.close();
				
				SAMRecord samRecord;
				//get relevant read sequences
				if(!prebuffer) {
					date = new Date();
					System.out.println(String.format("[%s]\tBuffering next read sequences.",date.toLocaleString()));
					readid2sequence.clear();
					recordIterator = bamReader.query(currentChrName, 0, 0, true);
					while(recordIterator.hasNext()) {
						samRecord = recordIterator.next();
						readid2sequence.put(samRecord.getReadName(),null);
					}
					recordIterator.close();
					sequenceReader = new BufferedReader(new FileReader(new File(fastaFilePath)));
					
					while((currentLine = sequenceReader.readLine()) != null) {
						if(currentLine.charAt(0) == '>' && readid2sequence.containsKey(currentLine.substring(1)))
							readid2sequence.put(currentLine.substring(1), sequenceReader.readLine());
					}
					sequenceReader.close();
				}
				
				//get and modifiy records from current chromosome
				date = new Date();
				System.out.println(String.format("[%s]\tModifying records from current chromosome.",date.toLocaleString()));
				recordIterator = bamReader.query(currentChrName, 0, 0, true);
				while(recordIterator.hasNext()) {
					samRecord = recordIterator.next();
					if(readid2sequence.get(samRecord.getReadName()) == null) {
						System.out.println(String.format("[%s]\tWarning: Cannot find sequence for read: %s.",date.toLocaleString(),samRecord.getReadName()));
						continue;
					}
					
					readSequence.setLength(0);
					readSequence.append(readid2sequence.get(samRecord.getReadName()).toUpperCase());
					modifyCigar(samRecord,chrSequence,readSequence,checkModification);
					bamWriter.addAlignment(samRecord);
				}
				recordIterator.close();
				System.out.println();
			}
			
			bamReader.close();
			bamWriter.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void modifyCigar(SAMRecord samRecord,StringBuffer chrSequence,StringBuffer readSequence,boolean checkModification) {
		Cigar cigar = samRecord.getCigar();
		Cigar modifiedCigar = new Cigar();
		
		//sam records are 1 based.
		int currentGenomePosition = samRecord.getAlignmentStart() - 1;
		int currentReadPosition = 0;
		
		if(samRecord.getReadNegativeStrandFlag())
			buildReverseComplement(readSequence);
		
		CigarElement currentElement;
		CigarOperator currentOperator;
		int currentMatchCount;
		int currentMismatchCount;
		for(int i = 0; i < cigar.numCigarElements(); i++) {
			currentElement = cigar.getCigarElement(i);
			currentOperator = currentElement.getOperator();
			if(currentOperator.name().equals("M")) {
				//check sequence here and replace M by a combination of '=' and 'X'
				currentMatchCount = 0;
				currentMismatchCount = 0;
			
				for(int j = 0; j < currentElement.getLength(); j++) {
					//found a match
					if(chrSequence.charAt(currentGenomePosition + j) == readSequence.charAt(currentReadPosition + j)) {
						//previous position was a mismatch
						if(currentMismatchCount != 0) {
							modifiedCigar.add(new CigarElement(currentMismatchCount, CigarOperator.X));
							currentMismatchCount = 0;
						}
						currentMatchCount++;
					}
					
					//found a mismatch
					else {
						if(currentMatchCount != 0) {
							modifiedCigar.add(new CigarElement(currentMatchCount,CigarOperator.EQ));
							currentMatchCount = 0;
						}
						currentMismatchCount++;
					}
				}
				//add last part
				if(currentMismatchCount != 0)
					modifiedCigar.add(new CigarElement(currentMismatchCount, CigarOperator.X));
				else if(currentMatchCount != 0)
					modifiedCigar.add(new CigarElement(currentMatchCount,CigarOperator.EQ));
			}
			else
				modifiedCigar.add(currentElement);
			
			//update read and genome position
			if(currentOperator.name().equals("M") || currentOperator.name().equals("I") || currentOperator.name().equals("X") || currentOperator.name().equals("=") || currentOperator.name().equals("S")) {
				currentGenomePosition += currentElement.getLength();
				currentReadPosition += currentElement.getLength();
			}
			if(currentOperator.name().equals("N")) {
				currentGenomePosition += currentElement.getLength();
			}
		}
		if(checkModification) checkModification(samRecord,modifiedCigar);
		samRecord.setCigar(modifiedCigar);
	}
	
	
	
	
	
	
	private void buildReverseComplement(StringBuffer sequence) {
		sequence.reverse();
		for(int i = 0; i < sequence.length(); i++) {
			if(sequence.charAt(i) == 'A' || sequence.charAt(i) == 'a')
				sequence.setCharAt(i, 'T');
			else if(sequence.charAt(i) == 'T' || sequence.charAt(i) == 't')
				sequence.setCharAt(i, 'A');
			else if(sequence.charAt(i) == 'G' || sequence.charAt(i) == 'g')
				sequence.setCharAt(i, 'C');
			else if(sequence.charAt(i) == 'C' || sequence.charAt(i) == 'c')
				sequence.setCharAt(i, 'G');
		}
	}
	
	
	private void checkModification(SAMRecord record, Cigar newCigar) {
		Date date;
		Cigar oldCigar = record.getCigar();
		int readLength = oldCigar.getReadLength();
		if(readLength != newCigar.getReadLength()) {
			date = new Date();
			System.err.println(String.format("[%s]\tWARNING: Cigar modification changed read length: %s -> %s",date.toLocaleString(),readLength,newCigar.getReadLength()));
			System.err.println(String.format("[%s]\tSAMLINE: %s",date.toLocaleString(),record.getSAMString()));
			System.err.println(String.format("[%s]\tNEWCIGAR: %s",date.toLocaleString(),newCigar.toString()));
		}
		
		try {
			int mismatches = record.getIntegerAttribute("NM");
			int newMismatches = 0;
			for(CigarElement element : newCigar.getCigarElements()) {
				if(element.getOperator().name() == "X")
					newMismatches += element.getLength();
			}
			if(mismatches != newMismatches) {
				date = new Date();
				System.err.println(String.format("[%s]\tWARNING: Cigar modification changed mismatch count: %s -> %s",date.toLocaleString(),mismatches,newMismatches));
				System.err.println(String.format("[%s]\tSAMLINE: %s",date.toLocaleString(),record.getSAMString()));
				System.err.println(String.format("[%s]\tNEWCIGAR: %s",date.toLocaleString(),newCigar.toString()));
			}
		}
		catch(Exception e) {
			date = new Date();
			System.err.println(String.format("[%s]\tWARNING: Sam file does not contain a field with NM tag. Can't check cigar modification.",date.toLocaleString()));
		}
	}
	
}
