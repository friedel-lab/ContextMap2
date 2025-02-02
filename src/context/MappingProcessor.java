package context;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NavigableMap;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Pattern;

import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;

import augmentedTree.IntervalTree;
import tools.BufferedRandomAccessFile;
import tools.String2Bitset;
import main.Context;
import main.InitialRead;
import main.InitialReadLocation;
import main.MutableInt;
import main.Pair;
import main.Read;
import main.ReadLocation;
import main.ReadPair;
import main.SparseRead;
import main.MultiSparseReadLocation;
import main.SparseReadLocation;
import main.Triplet;
import main.UnsynchronizedBufferedWriter;

public class MappingProcessor {
	
	
	//private HashMap<String,File> chr2file;
	private HashMap<String,File> chr2rmap;
	private ConcurrentMap<String,long[]> read2sequence;
	private String2Bitset string2bitset;
	private HashSet<String> spliceSites;
	private HashMap<String,TreeMap<Integer,HashSet<Integer>>> chr2annotatedSpliceSites;
	//in case we have given strand specific reads, chr2annotatedSpliceSites will be empty. instead
	//we split annotated splice sites by strand information
	private HashMap<String,TreeMap<Integer,HashSet<Integer>>> chr2annotatedForwardSpliceSites;
	private HashMap<String,TreeMap<Integer,HashSet<Integer>>> chr2annotatedReverseSpliceSites;
	
	private String readsPath;
	private String readFormat;
	private int readLength;
	private int maxMissmatches;
	private int maxMissmatchDifference;
	private int maxIntronLength;
	private int maxDelSize;
	private int seedLength;
	private int seedMissmatches;
	private int minPolyALength;
	private int minPolyAReadCount;
	private int maxConsideredClippingLength;
	private double upperCutoff;
	private double lowerCutoff;
	
	
	private StringBuilder currentChrSequence;
	
	
	private String multiMappingFilePath;
	private String outputPath;
	
	private boolean preferExtensionsWithKnownSpliceSignal;
	private boolean skipDenovoJunctions;
	private boolean skipNonCanonicalJunctions;
	private boolean clipping;
	private boolean polyA;
	private boolean strandedPolyA;
	
	
	
	public MappingProcessor() {
		this.string2bitset = new String2Bitset();
	}
	
	
	public MappingProcessor(String multiMappingFilePath, String outputPath, int maxMissmatchDifference,boolean preferExtensionsWithKnownSpliceSignal,boolean skipDenovoJunctions, boolean skipNonCanonicalJunctions) {
		this.multiMappingFilePath = multiMappingFilePath;
		this.outputPath = outputPath;
		this.maxMissmatchDifference = maxMissmatchDifference;
		this.preferExtensionsWithKnownSpliceSignal = preferExtensionsWithKnownSpliceSignal;
		this.skipDenovoJunctions = skipDenovoJunctions;
		this.skipNonCanonicalJunctions = skipNonCanonicalJunctions;
		this.string2bitset = new String2Bitset();
		
	}
	

	
	public MappingProcessor(String rmapFilesDir, String genomeDir, String annotationFilePath, String gtfFilePath, String readsPath, String readFormat, int readLength, int maxMissmatches, int maxMissmatchDifference, int maxIntronLength, int maxDelSize, int seedLength,int seedMissmatches, boolean preferExtensionsWithKnownSpliceSignal,boolean skipDenovoJunctions, boolean skipNonCanonicalJunctions, boolean strandSpecific, boolean clipping, boolean polyA, boolean strandedPolyA, int minPolyALength, int minPolyAReadCount, double upperCutoff, double lowerCutoff, int maxConsideredClippingLength) {
		this.readsPath = readsPath;
		this.readLength = readLength;
		this.maxMissmatches = maxMissmatches;
		this.maxMissmatchDifference = maxMissmatchDifference;
		this.maxIntronLength = maxIntronLength;
		this.maxDelSize = maxDelSize;
		this.seedLength = seedLength;
		this.seedMissmatches = seedMissmatches;
		this.preferExtensionsWithKnownSpliceSignal = preferExtensionsWithKnownSpliceSignal;
		this.skipDenovoJunctions = skipDenovoJunctions;
		this.skipNonCanonicalJunctions = skipNonCanonicalJunctions;
		this.clipping = clipping;
		this.polyA = polyA;
		this.strandedPolyA = strandedPolyA;
		this.minPolyALength = minPolyALength;
		this.minPolyAReadCount = minPolyAReadCount;
		this.upperCutoff = upperCutoff;
		this.lowerCutoff = lowerCutoff;
		this.maxConsideredClippingLength = maxConsideredClippingLength;
		this.string2bitset = new String2Bitset();
		if(readFormat == null)
			readFormat = "fasta";
				
		else
			this.readFormat = readFormat;
		
		try {
			
			//init the readers
			this.chr2rmap = new HashMap<String,File>();
			File[] rmapFiles = new File(rmapFilesDir).listFiles();
			String chr;
			for(File rmapFile : rmapFiles) {
				chr = rmapFile.getName().substring(0,rmapFile.getName().lastIndexOf("."));
				this.chr2rmap.put(chr, rmapFile);
			}
			
			//init the characteristic splice site patterns
			this.spliceSites = new HashSet<String>();
			//String[] characterizedSpliceSites = new String[]{"GTAG","GCAG","ATAC"};
			String[] characterizedSpliceSites = new String[]{"GTAG"};
			
			for(String spliceSite : characterizedSpliceSites) {
				this.spliceSites.add(spliceSite);
			}
			
			//in case we have given an annotation we hash known splice sites here
			if(annotationFilePath != null) {
				ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>> annotatedSpliceSites = parseAnnotatedSplitPositions(annotationFilePath,strandSpecific); 
				if(!strandSpecific)
					this.chr2annotatedSpliceSites = annotatedSpliceSites.get(0);
				else {
					this.chr2annotatedForwardSpliceSites = annotatedSpliceSites.get(0);
					this.chr2annotatedReverseSpliceSites = annotatedSpliceSites.get(1);
				}
			}
			else {
				this.chr2annotatedSpliceSites = null;
				this.chr2annotatedForwardSpliceSites = null;
				this.chr2annotatedReverseSpliceSites = null;
			}
			

			if(gtfFilePath != null) {
				ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>> annotatedSpliceSites = parseAnnotatedSplitPositionsFromGtf(gtfFilePath, strandSpecific);
				if(!strandSpecific)
					this.chr2annotatedSpliceSites = annotatedSpliceSites.get(0);
				else {
					this.chr2annotatedForwardSpliceSites = annotatedSpliceSites.get(0);
					this.chr2annotatedReverseSpliceSites = annotatedSpliceSites.get(1);
				}
			}
			
			this.currentChrSequence = new StringBuilder();
			DB db = DBMaker.newMemoryDB().transactionDisable().make();
			this.read2sequence = db.getHashMap("read2sequence");
			//this.read2sequence = new ConcurrentHashMap<String,String>();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void setCurrentChromosome(StringBuilder sb) {
		this.currentChrSequence = sb;
	}
	
	public StringBuilder getCurrentChromosome() {
		return this.currentChrSequence;
	}
	
	public void setCurrentReadSequences(ConcurrentMap<String,long[]> read2sequence) {
		this.read2sequence = read2sequence;
	}
	
	public ConcurrentMap<String,long[]> getCurrentReadSequence() {
		return this.read2sequence;
	}
	
	
	public void emptyReadSequences() {
		this.read2sequence.clear();
	}
	
	public void emptyChrSequence() {
		this.currentChrSequence.setLength(0);
	}
	
	public void parseMapping(Context context,BufferedRandomAccessFile rmapReader, boolean pairedEnd,StringBuilder localReference, int contextOffset, String alignmentsFilePath) {
		try {
			
			//first check if the actual chromosome sequence is set/available
			if(this.currentChrSequence.length() == 0) {
				System.err.println("Warning: Could not find sequence for chr:\t" + context.getChr());
				return;
			}
			
			if(alignmentsFilePath == null)
				parseReadsOfContext(context,localReference, contextOffset,rmapReader,pairedEnd);
			
			else
				parseReadsOfLargeContext(context,localReference,contextOffset,alignmentsFilePath);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public void extendMapping(Context context, StringBuilder localReference, int contextOffset, String multiMappingOutputPath, String mscsOutputFilePath, int minStartPositionOverlaps, String alignmentsFilePath) {
		try {
			//first check if the actual chromosome sequence is set/available
			if(this.currentChrSequence.length() == 0) {
				System.err.println("Warning: Could not find sequence for chr:\t" + context.getChr());
				return;
			}
			
			//generate a location pool which is used by the split extension and full extension method
			//this way we do not need to generate millions of ReadLocation Objects per Thread....
			ArrayList<ReadLocation> readLocationPool = new ArrayList<ReadLocation>(100);
			
			//extend the split positions
			UnsynchronizedBufferedWriter contextWriter = new UnsynchronizedBufferedWriter(new FileWriter(new File(multiMappingOutputPath)),1024 * 1024);
			ArrayList<Object> hashedSplitsAndFilePointer = determineSplitPositionsAndPrintFullReads(context,localReference,contextOffset,this.readsPath, this.readFormat,this.maxMissmatches,minStartPositionOverlaps,readLocationPool,contextWriter, alignmentsFilePath);
			TreeMap<Integer,HashSet<Integer>> hashedSplits = (TreeMap<Integer,HashSet<Integer>>)hashedSplitsAndFilePointer.get(0);
			long filePointer = (Long)hashedSplitsAndFilePointer.get(1);
			HashSet<String> insertionKeys = (HashSet<String>)hashedSplitsAndFilePointer.get(2);
			addAnnotatedSpliceSitesToHash(context,hashedSplits);
			
			
			
			//check if full and partial reads overlap with split positions and add these possible alignments
			//to the initial mapping
			extendPartialReads(context,localReference,contextOffset,hashedSplits,insertionKeys,filePointer,readLocationPool,contextWriter,alignmentsFilePath);
			contextWriter.close();
			hashedSplits.clear();
		}
		catch(Exception e) {
			e.printStackTrace();
			System.err.println("Error in local context: " + context.getId());
		}
	}
		
		
	private class ReadSequenceComparator implements Comparator<InitialRead> {
	
		private HashMap<String,String> read2sequence;
		
		public ReadSequenceComparator(ConcurrentHashMap<String,String> read2sequence,Context context) {

			this.read2sequence = new HashMap<String,String>();
			for(InitialRead read : context.getReads()) {
				this.read2sequence.put(read.getId(), read2sequence.get(read.getId()));
			}
		}
		
		public int compare(InitialRead l1, InitialRead l2) {
			String sequenceA = this.read2sequence.get(l1.getId());
			String sequenceB = this.read2sequence.get(l2.getId());
			return (sequenceA.compareTo(sequenceB));
		}
	}

	private void parseReadsOfContext(Context context,StringBuilder localReference, int contextOffset, BufferedRandomAccessFile rmapReader, boolean pairedEnd) {
		try {
			BufferedRandomAccessFile br;
			if(rmapReader == null) 
				br = new BufferedRandomAccessFile(this.chr2rmap.get(context.getChr()), "r",100000);
			else
				br = rmapReader;
				
			br.seek(context.getPointerToFirstRead());
			String currentLine;
			StringTokenizer st;
			String readId;
			String parentReadId;
			String readSequence;
			String mateSequence;
			String mateId = null;
			char mateInfo;
			char mappingType;
			String chr = context.getChr();
			int start;
			int end;
			int readLength = this.readLength;
			int mappingCount;
			int validPairsCount;
			String endAsString;
			String strand;
			char strandAsChar;
			
			int missmatches;
			boolean isMsc;
			
			String2Bitset s2b = new String2Bitset();
			
			DB db = DBMaker.newMemoryDB().transactionDisable().make();
			HTreeMap<String,Integer> read2index = db.createHashMap("read2index").keySerializer(Serializer.STRING).valueSerializer(Serializer.INTEGER).makeOrGet();
			HTreeMap<String,Integer> seq2readIndex = db.createHashMap("seq2readIndex").keySerializer(Serializer.STRING).valueSerializer(Serializer.INTEGER).makeOrGet();
			HTreeMap<String,Integer> mscSeq2readIndex = db.createHashMap("mscSeq2readIndex").keySerializer(Serializer.STRING).valueSerializer(Serializer.INTEGER).makeOrGet();
			HashSet<String> alreadyParsedSplits = new HashSet<String>();
			StringBuilder splitKeyBuilder = new StringBuilder();
			String splitKey;
			int currentIndex = 0;
			int mscIndex;
			InitialRead duplicatedRead;
			
			
			while((currentLine = br.getNextLine()) != null) {
				
				if(alreadyParsedSplits.size() >= 1000000)
					alreadyParsedSplits.clear();
				
				st = new StringTokenizer(currentLine,"\t");
				readId = st.nextToken();
				
				mappingType = st.nextToken().toCharArray()[0];
				st.nextToken();
				start = Integer.valueOf(st.nextToken());
				endAsString = st.nextToken();
				strand = st.nextToken();
				strandAsChar = strand.charAt(0);
				missmatches = Integer.valueOf(st.nextToken());
				readLength = Integer.valueOf(st.nextToken());
				mappingCount = Integer.valueOf(st.nextToken());
				validPairsCount = Integer.valueOf(st.nextToken());
				
				
				parentReadId = readId;
				isMsc = false;
				duplicatedRead = null;
				mscIndex = readId.indexOf("::MSC");
				readSequence = null;
				if(!read2index.containsKey(readId) || mscIndex != -1) {
					readSequence = s2b.decompress(this.read2sequence.get(readId));
				}
				
				
				mateSequence = null;
				if(pairedEnd) {
					mateInfo = readId.charAt(readId.length() - 1);
					if(mateInfo == '1')
						mateId = readId.substring(0,readId.lastIndexOf('/')) + "/2";
					else
						mateId = readId.substring(0,readId.lastIndexOf('/')) + "/1";
					
					if(!isMsc && this.read2sequence.containsKey(mateId)) {
						mateSequence = s2b.decompress(this.read2sequence.get(mateId));
					}
					
				}
				
				
				if(mscIndex != -1) {
					parentReadId = readId.substring(0,mscIndex);
					isMsc = true;
					if(mscSeq2readIndex.containsKey(readSequence))
						duplicatedRead = context.getRead(mscSeq2readIndex.get(readSequence));
					
				}
				
				splitKey = null;
				if(mappingType == 'S') {
					end = Integer.valueOf(endAsString);
					splitKeyBuilder.setLength(0);
					splitKey = splitKeyBuilder.append(readId).append('_').append(start).append('_').append(end).append('_').append(strand).toString();
				}
				else {
					end = start + readLength - 1;
				}
				

				if((!context.isStrandSpecific() || context.getStrand().equals(strand)) && start > 0 && end - contextOffset <= localReference.length()) {
					
					if(read2index.containsKey(readId)) {
						if(mappingType != 'S' || !alreadyParsedSplits.contains(splitKey)) {
							context.getRead(read2index.get(readId)).addLocation(strandAsChar, mappingType, start, end, missmatches);
						}
					}
					
					//no msc, not pairedend or no valid pair found or mate not contained in context
					else if((!pairedEnd || validPairsCount == 0) && !isMsc && seq2readIndex.containsKey(readSequence)) {
						context.getRead(seq2readIndex.get(readSequence)).addDuplicate(readId);
					}
					
					//no msc, pairedend and sequences of both mates equal
					else if(pairedEnd && !isMsc && mateSequence != null && seq2readIndex.containsKey(readSequence + '_' + mateSequence)) {
						context.getRead(seq2readIndex.get(readSequence + '_' + mateSequence)).addDuplicate(readId);
						
					}
					
					//msc, not pairedend or no valid pair found
					else if((!pairedEnd || validPairsCount == 0) && isMsc && mscSeq2readIndex.containsKey(readSequence) && s2b.decompress(this.read2sequence.get(parentReadId)).equals(s2b.decompress(this.read2sequence.get(duplicatedRead.getId().substring(0,duplicatedRead.getId().indexOf("::MSC")))))) {
						duplicatedRead.addDuplicate(readId);
						
					}
					
					else {
						InitialRead read = new InitialRead(readId, strandAsChar, mappingType, start, end, missmatches,mappingCount,validPairsCount);
						context.addRead(read);
						read2index.put(readId, currentIndex);
						
						
						
						if(!isMsc) {
							if(!pairedEnd || validPairsCount == 0)
								seq2readIndex.put(readSequence, currentIndex);
							else if(pairedEnd && mateSequence != null)
								seq2readIndex.put(readSequence + '_' + mateSequence,currentIndex);
								
						}
						else {
							if(!pairedEnd || validPairsCount == 0) {
								mscSeq2readIndex.put(readSequence,currentIndex);
							}
						}
						currentIndex++;
					}
					
					if(mappingType == 'S')
						alreadyParsedSplits.add(splitKey);
					
				}
				
				if(br.getFilePointer() == context.getPointerToLastRead()) {
					break;
				}
				
				}
				if(rmapReader == null) 
					br.close();
				
				
				db.close();
			}
		
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	
	private void parseReadsOfLargeContext(Context context,StringBuilder localReference, int contextOffset, String alignmentsFilePath) {
		try {
			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(alignmentsFilePath), "r",10240);
				
			String currentLine;
			StringTokenizer st;
			String readId;
			String prevReadId = "";
			char mappingType;
			int start;
			int end;
			int readLength = this.readLength;
			int mappingCount;
			int validPairsCount;
			String endAsString;
			String strand;
			char strandAsChar;
			int missmatches;
			HashSet<String> alreadyParsedSplits = new HashSet<String>();
			StringBuilder splitKeyBuilder = new StringBuilder();
			String splitKey;
			long prevPointer= braf.getFilePointer();
			while((currentLine = braf.getNextLine()) != null) {
				st = new StringTokenizer(currentLine,"\t");
				readId = st.nextToken();
				mappingType = st.nextToken().toCharArray()[0];
				st.nextToken();
				start = Integer.valueOf(st.nextToken());
				endAsString = st.nextToken();
				strand = st.nextToken();
				strandAsChar = strand.charAt(0);
				missmatches = Integer.valueOf(st.nextToken());
				readLength = Integer.valueOf(st.nextToken());
				mappingCount = Integer.valueOf(st.nextToken());
				validPairsCount = Integer.valueOf(st.nextToken());
				
				splitKey = null;
				if(mappingType == 'S') {
					end = Integer.valueOf(endAsString);
					splitKeyBuilder.setLength(0);
					splitKey = splitKeyBuilder.append(readId).append('_').append(start).append('_').append(end).append('_').append(strand).toString();
				}
				else {
					end = start + readLength - 1;
				}
				

				if((!context.isStrandSpecific() || context.getStrand().equals(strand)) && start > 0 && end - contextOffset <= localReference.length()) {
					
					if(readId.equals(prevReadId)) {
						//do nothing...
					}
					
					else {
						InitialRead read = new InitialRead(readId, prevPointer,mappingCount,validPairsCount);
						context.addRead(read);
						prevReadId = readId;
					}
				}
				
				prevPointer= braf.getFilePointer();
			}
			
			braf.close();
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	private void parseReadLocations(Context context,InitialRead read, int contextOffset, int localReferenceLength, BufferedRandomAccessFile braf) throws Exception {
		String currentLine;
		String readId;
		char mappingType;
		int start;
		String endAsString;
		int end;
		char strand;
		int mismatches;
		int readLength;
		StringTokenizer st;
		
		HashSet<String> alreadyParsedSplits = new HashSet<String>();
		StringBuilder splitKeyBuilder = new StringBuilder();
		String splitKey;
		
		ArrayList<InitialReadLocation> locations = new ArrayList<InitialReadLocation>();
		braf.seek(read.getLocationStartPointer());
		while((currentLine = braf.getNextLine()) != null) {
			st = new StringTokenizer(currentLine,"\t");
			readId = st.nextToken();
			if(!readId.equals(read.getId()))
				break;
			
			mappingType = st.nextToken().toCharArray()[0];
			st.nextToken();
			start = Integer.valueOf(st.nextToken());
			endAsString = st.nextToken();
			strand = st.nextToken().toCharArray()[0];
			mismatches = Integer.valueOf(st.nextToken());
			readLength = Integer.valueOf(st.nextToken());
			
			splitKey = null;
			if(mappingType == 'S') {
				end = Integer.valueOf(endAsString);
				splitKeyBuilder.setLength(0);
				splitKey = splitKeyBuilder.append(readId).append('_').append(start).append('_').append(end).append('_').append(strand).toString();
			}
			else {
				end = start + readLength - 1;
			}
			
			if((!context.isStrandSpecific() || context.getStrand().equals(strand)) && start > 0 && end - contextOffset <= localReferenceLength) {
				if(mappingType != 'S' || !alreadyParsedSplits.contains(splitKey)) {
					locations.add(new InitialReadLocation(strand,mappingType, start, end, mismatches));
					if(mappingType == 'S')
						alreadyParsedSplits.add(splitKey);
				}
			}
		}
		
		read.setLocations(locations);
	}
	
	private ArrayList<Object> determineSplitPositionsAndPrintFullReads(Context context,StringBuilder localReference,int contextOffset, String readsPath, String readFormat, int maxMissmatches,int minStartPositionOverlaps, ArrayList<ReadLocation> readLocationPool, UnsynchronizedBufferedWriter pw, String alignmentsFilePath) throws Exception {
		ArrayList<Object> returnValues = new ArrayList<Object>();
		HashMap<String,Split> key2split = new HashMap<String,Split>();
		String splitKey;
		TreeMap<Integer,HashSet<Integer>> hashedSplits = new TreeMap<Integer,HashSet<Integer>>();
		HashSet<String> insertionKeys = new HashSet<String>();
		HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentAnnotatedSpliceSiteMap = null;
		if(!context.isStrandSpecific())
			currentAnnotatedSpliceSiteMap = this.chr2annotatedSpliceSites;
		if(context.isStrandSpecific() && context.getStrand().equals("+"))
			currentAnnotatedSpliceSiteMap = this.chr2annotatedForwardSpliceSites;
		else if(context.isStrandSpecific() && context.getStrand().equals("-"))
			currentAnnotatedSpliceSiteMap = this.chr2annotatedReverseSpliceSites;
		
		int startA = -1;
		int endA = -1;
		int startB = -1;
		int endB = -1;
		int mismatches = -1;
		int minMismatches = -1;
		int minMismatchesOfFullAlignment = -1;
		int insertionSize = 0;
		long filePointer = 0;
		boolean hasInsertion = false;
		char strandOfSpliceSignal = '0';
		boolean isIndel;
		boolean foundLocationWithKnownSpliceSignal = false;
		boolean foundLocationOverlappingKnownJunction = false;
		boolean bestMappingWithoutSignal = false;
		boolean foundFullReadAlignment = false;
		boolean foundPartialAlignment = false;
		ArrayList<ReadLocation> locationsToAdd = new ArrayList<ReadLocation>();
		HashSet<ReadLocation> locationsToRemove = new HashSet<ReadLocation>();
		String readSequence;
		String prevReadSequence = " ";
		StringBuilder readSequenceBuffer = new StringBuilder();
		StringBuilder donorSite = new StringBuilder();
		StringBuilder acceptorSite = new StringBuilder();
		StringBuilder tmpLineBuilder = new StringBuilder();
		int[] mismatchesSplitPartA;
		int[] mismatchesSplitPartB;
		ReadLocation tmpLocation;
		String tmpLine;
		HashMap<Integer,MutableInt> splitSize2frequency = new HashMap<Integer,MutableInt>();
		int splitSize;
		
		//here we know that the context is not in memory yet. we just parse the alignments of the current read.
		BufferedRandomAccessFile braf = null;
		if(alignmentsFilePath != null) {
			braf = new BufferedRandomAccessFile(new File(alignmentsFilePath),"r",1024);
		}
		
		for(InitialRead read : context.getReads()) {
			readSequence = this.string2bitset.decompress(this.read2sequence.get(read.getId()));
			
			if(readSequence == null) {
				System.err.println(read.getId());
				System.exit(1);
			}
			
			
			if(read.getLocations() == null && braf != null) {
				parseReadLocations(context,read,contextOffset,localReference.length(),braf);
			}
			
			mismatchesSplitPartA = new int[readSequence.length()];
			mismatchesSplitPartB = new int[readSequence.length()];
		
			prevReadSequence = readSequence;
			locationsToAdd.clear();
			minMismatches = this.maxMissmatches + 1;
			minMismatchesOfFullAlignment = this.maxMissmatches + 1;
			bestMappingWithoutSignal = true;
			foundFullReadAlignment = false;
			foundPartialAlignment = false;
			for(InitialReadLocation readLocation : read.getLocations()) {
				if(readLocation.getMappingType() == 'S') {
					
					startA = readLocation.getStartA();
					endB = readLocation.getEndA();

					getMismatchesForSplitPartA(mismatchesSplitPartA,readSequence,localReference,contextOffset,readSequenceBuffer,startA,startA + readSequence.length() - 1,readLocation.getStrand());
					getMismatchesForSplitPartB(mismatchesSplitPartB,readSequence,localReference,contextOffset,readSequenceBuffer,endB - readSequence.length() + 1,endB,readLocation.getStrand());
					//j defines the length of the first split part
					for(int j = 1; j < readSequence.length(); j++) {
						endA = startA + j - 1;
						startB = endB - (readSequence.length() - j) + 1;
						
						//splits with negative "intron length" -> insertion
						hasInsertion = false;
						if(startA + readSequence.length() - 1 > endB) {
							insertionSize = (startA + readSequence.length() - 1) - endB;
							
							if(readSequence.length() - j - 1 - insertionSize < 0 || (j + insertionSize == readSequence.length()))
								break;
							
							mismatches = mismatchesSplitPartA[j-1] + mismatchesSplitPartB[readSequence.length() - j - 1 - insertionSize];
							hasInsertion = true;
						}
						
						else {
							mismatches = mismatchesSplitPartA[j-1] + mismatchesSplitPartB[readSequence.length() - j - 1];
						}
						
						
						if(mismatches <= maxMissmatches) {

							if(locationsToAdd.size() < readLocationPool.size()) {
								tmpLocation = readLocationPool.get(locationsToAdd.size());
								tmpLocation.updateLocation(context.getChr(), readLocation.getStrand(),'S',startA,endA,startB,endB,mismatches);
							}
							
							else {
								tmpLocation = new ReadLocation(context.getChr(), readLocation.getStrand(),'S',startA,endA,startB,endB,mismatches);
								readLocationPool.add(tmpLocation);
							}
							
							if(!hasInsertion && (startB - endA - 1) > this.maxDelSize) {
								//check if the actual split has a known splice signal
								strandOfSpliceSignal = hasSpliceSignal(localReference,contextOffset, donorSite, acceptorSite,endA,startB,tmpLocation.getStrand(),context.isStrandSpecific());
								tmpLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
								
								
								//in case an annotation is given check if the actual split overlaps with an annotated junction
								if(currentAnnotatedSpliceSiteMap != null && currentAnnotatedSpliceSiteMap.containsKey(tmpLocation.getChr())) {
									if((currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getEndA()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getEndA()).contains(tmpLocation.getStartB())) ||
									   (currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getStartB()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getStartB()).contains(tmpLocation.getEndA()))){
										tmpLocation.setOverlapsKnownJunction(true);
										
									}
								}
							}
							
							if(mismatches < minMismatches) {
								minMismatches = mismatches;
								if(tmpLocation.hasSpliceSignal() || tmpLocation.overlapsKnownJunction()) {
									bestMappingWithoutSignal = false;
								}
								else
									bestMappingWithoutSignal = true;
							}
							
							else if(mismatches == minMismatches) {
								if(tmpLocation.hasSpliceSignal() || tmpLocation.overlapsKnownJunction())
									bestMappingWithoutSignal = false;
							}
							
							
							locationsToAdd.add(tmpLocation);
							
						}
						else {
							j += (mismatches - this.maxMissmatches) - 1;
						}
					}
				
					//check last location
					if(mismatches <= maxMissmatches) {
			
						if(locationsToAdd.size() < readLocationPool.size()) {
							tmpLocation = readLocationPool.get(locationsToAdd.size());
							tmpLocation.updateLocation(context.getChr(), readLocation.getStrand(),'S',startA,endA,startB,endB,mismatches);
						}
						
						else {
							tmpLocation = new ReadLocation(context.getChr(), readLocation.getStrand(),'S',startA,endA,startB,endB,mismatches);
							readLocationPool.add(tmpLocation);
						}
						
						
						if(!hasInsertion && (startB - endA - 1) > this.maxDelSize) {
							strandOfSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA,startB,tmpLocation.getStrand(),context.isStrandSpecific());
							tmpLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
							
							if(currentAnnotatedSpliceSiteMap != null && currentAnnotatedSpliceSiteMap.containsKey(tmpLocation.getChr())) {
								if((currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getEndA()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getEndA()).contains(tmpLocation.getStartB())) ||
								   (currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getStartB()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getStartB()).contains(tmpLocation.getEndA()))){
									tmpLocation.setOverlapsKnownJunction(true);
								}
							}
						}
						
						if(mismatches < minMismatches) {
							minMismatches = mismatches;
							if(tmpLocation.hasSpliceSignal() || tmpLocation.overlapsKnownJunction()) {
								bestMappingWithoutSignal = false;
							}
							else
								bestMappingWithoutSignal = true;
						}
						
						else if(mismatches == minMismatches) {
							if(tmpLocation.hasSpliceSignal() || tmpLocation.overlapsKnownJunction())
								bestMappingWithoutSignal = false;
						}
						
						
						
						locationsToAdd.add(tmpLocation);
					}
				}
				else if(readLocation.getMappingType() == 'F') {
					foundFullReadAlignment = true;
					if(readLocation.getMissmatches() < minMismatchesOfFullAlignment)
						minMismatchesOfFullAlignment = readLocation.getMissmatches();
				}
				
				else if(readLocation.getMappingType() == 'P') {
					foundPartialAlignment = true;
				}
			}
			
			foundLocationWithKnownSpliceSignal = false;
			foundLocationOverlappingKnownJunction = false;
			for(ReadLocation rl : locationsToAdd) {
				if(rl.getMismatches() - minMismatches <= 1) {
					if(rl.hasSpliceSignal())
						foundLocationWithKnownSpliceSignal = true;
					if(rl.overlapsKnownJunction())
						foundLocationOverlappingKnownJunction = true;
				}
			}
				
			
			
			//first print the full reads of the current read
			filePointer += printFullReads(read,null,context.getId(),context.getChr(),pw);
			
			
			//now filter the split candidates
			locationsToRemove.clear();
			
			if(bestMappingWithoutSignal && (foundLocationWithKnownSpliceSignal || foundLocationOverlappingKnownJunction)) {
				
				int allowedMismatchDifferenceFromBestMapping;
				for(ReadLocation rl : locationsToAdd) {
					allowedMismatchDifferenceFromBestMapping = 0;
					if(bestMappingWithoutSignal && (rl.hasSpliceSignal() || rl.overlapsKnownJunction()))
						allowedMismatchDifferenceFromBestMapping = 1;
					
					if(rl.getMismatches() - minMismatches > allowedMismatchDifferenceFromBestMapping) {
						locationsToRemove.add(rl);
						continue;
					}
					if(foundLocationOverlappingKnownJunction) {
						if(!rl.overlapsKnownJunction()) {
							locationsToRemove.add(rl);
						}
						continue;
					}
					if(foundLocationWithKnownSpliceSignal) {
						if(!rl.hasSpliceSignal())
							locationsToRemove.add(rl);
						continue;
					}
				}
			}
			
			//here we know the best mapping comes from a known splice signal or we did not find any candidate with a signal. just remove those with too many mismatches
			else {
				for(ReadLocation rl : locationsToAdd) {
					if(rl.getMismatches() - minMismatches > this.maxMissmatchDifference) {
						locationsToRemove.add(rl);
					}
				}
			}
			
			
			//remove non canonical splits for reads with an full alignment with mismatch value <= split mismatch value or if the skip option is set, remove all those predictions
			for(ReadLocation rl : locationsToAdd) {
				if(!rl.overlapsKnownJunction() && !rl.hasSpliceSignal() && ((foundFullReadAlignment && minMismatchesOfFullAlignment <= rl.getMismatches()) || this.skipNonCanonicalJunctions)) {
					//check if we process an indel
					if(rl.getStartB() > rl.getEndA() && rl.getStartB() - rl.getEndA() > this.maxDelSize)
						locationsToRemove.add(rl);
				}
			}
			
			
			locationsToAdd.removeAll(locationsToRemove);
			
			
			
			
			if((foundFullReadAlignment && locationsToAdd.size() >= 2 && minMismatches >= minMismatchesOfFullAlignment) ||
			   (!foundFullReadAlignment && locationsToAdd.size() > 2)) {
				splitSize2frequency.clear();
				for(ReadLocation l : locationsToAdd) {
					splitSize = l.getEndA() - l.getStartA() + 1;
					if(splitSize2frequency.containsKey(splitSize))
						splitSize2frequency.get(splitSize).increment();
					else
						splitSize2frequency.put(splitSize, new MutableInt(1));
				}
				
				for(MutableInt freq : splitSize2frequency.values()) {
					if(freq.intValue() > locationsToAdd.size()/2) {
						locationsToAdd.clear();
						break;
					}
				}
			}
			
			
			
			
			//now print the split candidates
			for(ReadLocation locationToAdd : locationsToAdd) {
				if((currentAnnotatedSpliceSiteMap == null || !foundLocationOverlappingKnownJunction) && !this.skipDenovoJunctions) {
					if(!this.preferExtensionsWithKnownSpliceSignal || !foundLocationWithKnownSpliceSignal || locationToAdd.hasSpliceSignal()) {
						tmpLineBuilder.setLength(0);
						tmpLine = tmpLineBuilder.append(context.getId()).append("\t").append(read.getId()).append("\t").append(locationToAdd.getChr()).append("\t").append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
						pw.write(tmpLine);
						pw.newLine();
						filePointer += tmpLine.length() + 1;
						
						
						
						isIndel = (locationToAdd.getStartB() < locationToAdd.getEndA() || locationToAdd.getStartB() - locationToAdd.getEndA() <= this.maxDelSize);
						splitKey = locationToAdd.getEndA() + "_" + locationToAdd.getStartB();
						if(!key2split.containsKey(splitKey)) {
							key2split.put(splitKey, new Split(locationToAdd.getEndA(),locationToAdd.getStartB(),locationToAdd.hasSpliceSignal(),locationToAdd.overlapsKnownJunction(),isIndel));
						}
						key2split.get(splitKey).addStartPosition(locationToAdd.getStartA());
						
					}
				}
				
				else if(locationToAdd.overlapsKnownJunction()) {
					tmpLineBuilder.setLength(0);
					tmpLine = tmpLineBuilder.append(context.getId()).append("\t").append(read.getId()).append("\t").append(locationToAdd.getChr()).append("\t").append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
					pw.write(tmpLine);
					pw.newLine();
					filePointer += tmpLine.length() + 1;
					
					
					
					splitKey = locationToAdd.getEndA() + "_" + locationToAdd.getStartB();
					if(!key2split.containsKey(splitKey)) {
						key2split.put(splitKey, new Split(locationToAdd.getEndA(),locationToAdd.getStartB(),locationToAdd.hasSpliceSignal(),locationToAdd.overlapsKnownJunction(),false));
					}
					key2split.get(splitKey).addStartPosition(locationToAdd.getStartA());
				
				
				}
			}
			
			
			//now print the same for the duplicates
			if(read.getDuplicates() != null) {
				for(String id : read.getDuplicates()) {
					filePointer += printFullReads(read,id,context.getId(),context.getChr(),pw);
					for(ReadLocation locationToAdd : locationsToAdd) {
						if((currentAnnotatedSpliceSiteMap == null || !foundLocationOverlappingKnownJunction) && !this.skipDenovoJunctions) {
							if(!this.preferExtensionsWithKnownSpliceSignal || !foundLocationWithKnownSpliceSignal || locationToAdd.hasSpliceSignal()) {
								tmpLineBuilder.setLength(0);
								tmpLine = tmpLineBuilder.append(context.getId()).append("\t").append(id).append("\t").append(locationToAdd.getChr()).append("\t").append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
								pw.write(tmpLine);
								pw.newLine();
								filePointer += tmpLine.length() + 1;
							}
						}
						
						else if(locationToAdd.overlapsKnownJunction()) {
							tmpLineBuilder.setLength(0);
							tmpLine = tmpLineBuilder.append(context.getId()).append("\t").append(id).append("\t").append(locationToAdd.getChr()).append("\t").append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
							pw.write(tmpLine);
							pw.newLine();
							filePointer += tmpLine.length() + 1;
							
						}
					}
				}
			}
			
			
			//if there is no partial alignment available we can remove all locations
			if(!foundPartialAlignment) {
				read.getLocations().clear();
				if(read.getDuplicates() != null)
					read.getDuplicates().clear();
			}
			
			if(braf != null) {
				read.getLocations().clear();
				read.setLocations(null);
			}
		}
		
		
		
		//add splits to hash
		for(Split split : key2split.values()) {
			
			if(split.getDifferentStartPositions() < minStartPositionOverlaps)
				continue;
			
			
			if(hashedSplits.containsKey(split.getSplitPositionA()))
				hashedSplits.get(split.getSplitPositionA()).add(split.getSplitPositionB());
			else {
				HashSet<Integer> tmpSet = new HashSet<Integer>();
				tmpSet.add(split.getSplitPositionB());
				hashedSplits.put(split.getSplitPositionA(), tmpSet);
			}
			
			
			if(split.getSplitPositionA() < split.getSplitPositionB()) {
				if(hashedSplits.containsKey(split.getSplitPositionB())) {
					hashedSplits.get(split.getSplitPositionB()).add(split.getSplitPositionA());
				}
				else {
					HashSet<Integer> tmpSet = new HashSet<Integer>();
					tmpSet.add(split.getSplitPositionA());
					hashedSplits.put(split.getSplitPositionB(), tmpSet);
				}
			}
			
			else {
				//insertionKeys.add(String.format("%s_%s",split.getSplitPositionA(), split.getSplitPositionB()));
				insertionKeys.add(split.getSplitPositionA() + "_" + split.getSplitPositionB());
			}
			
		}
		
		
		if(braf != null) braf.close();
		
		returnValues.add(hashedSplits);
		returnValues.add(filePointer);
		returnValues.add(insertionKeys);
		return returnValues;
	}
	
	
	
	/**
	 * duplicateId will only be set for duplicated reads, otherwise set to 'null'
	 * @param read
	 * @param duplicateId
	 * @param contextId
	 * @param writer
	 * @return
	 * @throws Exception
	 */
	
	private long printFullReads(InitialRead read,String duplicateId, String contextId, String chr, UnsynchronizedBufferedWriter writer) throws Exception {
		
		long filePointer = 0;
		String currentLine;
		StringBuilder lineBuilder = new StringBuilder();
		String id = read.getId();
		if(duplicateId != null)
			id = duplicateId;
		for(InitialReadLocation location : read.getLocations()) {
			if(location.getMappingType() == 'F') {
				lineBuilder.setLength(0);
				currentLine = lineBuilder.append(contextId).append("\t").append(id).append("\t").append(chr).append("\t").append(location.getStrand()).append("\t").append(location.getStartA()).append("\t").append(location.getEndA()).append("\t").append(location.getStartB()).append("\t").append(location.getEndB()).append("\t").append(location.getMissmatches()).append("\t").append("0").append("\t").append(location.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
				//currentLine = lineBuilder.append(contextId).append("\t").append(id).append("\t").append(location.getChr()).append("\t").append(location.getStrand()).append("\t").append(location.getStartA()).append("\t").append(location.getEndA()).append("\t").append(location.getStartB()).append("\t").append(location.getEndB()).append("\t").append(location.getMissmatches()).append("\t").append("0").append("\t").append(location.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
				 
				writer.write(currentLine);
				writer.newLine();
				filePointer += currentLine.length() + 1;
			}
		}
		return filePointer;
	}
	
	
	private void getMismatchesForSplitPartA(int[] mismatches,String readSequence,StringBuilder localReference, int contextOffset, StringBuilder readSequenceBuffer, int start, int end, char strand) {
		int mismatchCount = 0;
		
		if(start - contextOffset < 0 || end - contextOffset + 1 > localReference.length()) {
			for(int i = 0; i < mismatches.length; i++) {
				mismatches[i] = this.maxMissmatches + 1;
			}
			return;
		}
		
		int referenceOffset = start - contextOffset;
		
		if(strand == '-') {
			readSequenceBuffer.setLength(0);
			readSequenceBuffer.append(readSequence);
			readSequenceBuffer.reverse();
			for(int i = 0; i < readSequenceBuffer.length(); i++) {
				readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
			}
			readSequence = readSequenceBuffer.toString();
		}
		
		int j = readSequence.length();
		for(int i = 0; i < readSequence.length(); i++) {
			if(readSequence.charAt(i) != localReference.charAt(referenceOffset + i))
				mismatchCount++;
			mismatches[i] = mismatchCount;
			
			if(mismatchCount > this.maxMissmatches) {
				j = i+1;
				break;
			}
		}
		
		for(int i = j; i < readSequence.length(); i++) {
			mismatches[i] = mismatchCount;
		}

	}
	
	/**
	 * 
	 * use strand '0' if the read sequence of an alignment to the negative strand is already reverse complemented. 
	 * 
	 * @param mismatches
	 * @param readSequence
	 * @param localReference
	 * @param contextOffset
	 * @param readSequenceBuffer
	 * @param start
	 * @param end
	 * @param strand
	 */
	
	private void getMismatchesForSplitPartB(int[] mismatches,String readSequence,StringBuilder localReference, int contextOffset, StringBuilder readSequenceBuffer, int start, int end, char strand) {
		int mismatchCount = 0;
		
		if(start - contextOffset <= 0 || end - contextOffset + 1 > localReference.length()) {
			for(int i = 0; i < mismatches.length; i++) {
				mismatches[i] = this.maxMissmatches + 1;
			}
			return;
		}
		
		int referenceOffset = start - contextOffset;
		
		if(strand == '-') {
			readSequenceBuffer.setLength(0);
			readSequenceBuffer.append(readSequence);
			readSequenceBuffer.reverse();
			for(int i = 0; i < readSequenceBuffer.length(); i++) {
				readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
			}
			readSequence = readSequenceBuffer.toString();
		}
		
		int j = -1;
		for(int i = readSequence.length() - 1; i >= 0; i--) {
			if(readSequence.charAt(i) != localReference.charAt(referenceOffset + i))
				mismatchCount++;
			mismatches[readSequence.length() -i -1] = mismatchCount;
			
			if(mismatchCount > this.maxMissmatches) {
				j = i - 1;
				break;
			}
		}
		
		for(int i = j; i >= 0; i--) {
			mismatches[readSequence.length() -i -1] = mismatchCount;
		}
	}
	
	
		
	
	private ArrayList<Triplet<Integer,Integer,Integer>> getBestLocalAlignmentIndices(String readSequence, StringBuilder localReference, int contextOffset, StringBuilder readSequenceBuffer, int start, int end, char strand, int matchScore, int mismatchPenalty, int maxMismatches, int minAlignmentLength, boolean clipAtStartOnly, boolean clipAtEndOnly) {
		ArrayList<Triplet<Integer,Integer,Integer>> results = new ArrayList<Triplet<Integer,Integer,Integer>>();
		
		if(start - contextOffset < 0 || (end - contextOffset + 1) > localReference.length())
			return results;
		
		int refOffset = start - contextOffset;
		//String referenceSequence = localReference.substring(start - contextOffset,end - contextOffset + 1).toUpperCase();
		
		if(strand == '-') {
			readSequenceBuffer.setLength(0);
			readSequenceBuffer.append(readSequence);
			readSequenceBuffer.reverse();
			for(int i = 0; i < readSequenceBuffer.length(); i++) {
				readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
			}
			readSequence = readSequenceBuffer.toString();
		}
		
		
		
		int currentScore = 0;
		int scoreMax = mismatchPenalty;
		int[] scores = new int[readSequence.length() + 1];
		scores[0] = 0;
		ArrayList<Integer> maxima = new ArrayList<Integer>();
		for(int i = 0; i < readSequence.length(); i++) {
			
			//update the score
			//if(readSequence.charAt(i) == referenceSequence.charAt(i)) {
			if(readSequence.charAt(i) == localReference.charAt(refOffset + i)) {
				currentScore += matchScore;
			}
			else {
				currentScore += mismatchPenalty;
			}
			
			//set the score
			if(currentScore < 0) {
				scores[i+1] = 0;
				currentScore = 0;
			}
			else
				scores[i+1] = currentScore;
			

			//update the max score
			if(currentScore > scoreMax) {
				scoreMax = currentScore;
				maxima.clear();
				maxima.add(i+1);
			}
			else if(currentScore == scoreMax) {
				maxima.add(i+1);
			}
			
		}
		
		//now go from every maxima to the next zero entry
		int alignmentLength;
		int mismatchCount;
		for(int maximumIndex : maxima) {
			for(int i = maximumIndex - 1; i>= 0; i--) {
				if(scores[i] == 0) {
					alignmentLength = maximumIndex - i;
					mismatchCount = ((alignmentLength * matchScore) - scores[maximumIndex])/(-1 * mismatchPenalty);
					if(alignmentLength >= minAlignmentLength && alignmentLength < readSequence.length() && mismatchCount <= maxMismatches) {
						
						if(clipAtStartOnly && (maximumIndex - 1) != readSequence.length() - 1)
							continue;
						
						if(clipAtEndOnly && i != 0)
							continue;
						
						Triplet<Integer,Integer,Integer> tmpResult = new Triplet<Integer,Integer,Integer>();
						tmpResult.setFirst(i);
						tmpResult.setSecond(maximumIndex - 1);
						tmpResult.setThird(mismatchCount);
						results.add(tmpResult);
					}
					
					break;
				}
			}
			
		}
		
		return results;
	}
	

	private char substitute(char n) {
		if(n == 'A' || n == 'a') return 'T';
		if(n == 'T' || n == 't') return 'A';
		if(n == 'C' || n == 'c') return 'G';
		if(n == 'G' || n == 'g') return 'C';
		else return 'N';
	}
	

	private char hasSpliceSignal(StringBuilder localReference, int contextOffset,StringBuilder donorSite, StringBuilder acceptorSite, int splitPointA, int splitPointB, char strand, boolean strandSpecific) {
		donorSite.setLength(0);
		acceptorSite.setLength(0);
		donorSite.append(localReference.substring((splitPointA + 1) - contextOffset, (splitPointA + 3) - contextOffset));
		acceptorSite.append(localReference.substring((splitPointB - 2) - contextOffset, splitPointB - contextOffset));
		
		//check forward splice sites
		if(!strandSpecific || strand == '+') {
		if(this.spliceSites.contains(donorSite.toString().toUpperCase() + acceptorSite.toString().toUpperCase()))
			return '+';
		}
		
		
		//check reverse splice sites
		donorSite.reverse();
		acceptorSite.reverse();
		for(int i = 0; i < acceptorSite.length(); i++) {
			acceptorSite.setCharAt(i, substitute(acceptorSite.charAt(i)));
			donorSite.setCharAt(i, substitute(donorSite.charAt(i)));
		}
		if(!strandSpecific || strand == '-') {
		if(this.spliceSites.contains(acceptorSite.toString().toUpperCase() + donorSite.toString().toUpperCase()))
			return '-';
		}
		
		return '0';
	}
	
			
	private void extendPartialReads(Context context,StringBuilder localReference,int contextOffset,TreeMap<Integer,HashSet<Integer>> hashedSplits, HashSet<String> insertionKeys, long filePointer, ArrayList<ReadLocation> readLocationPool, UnsynchronizedBufferedWriter pw, String alignmentsFilePath) {
		try {
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentAnnotatedSpliceSiteMap = null;
			if(!context.isStrandSpecific())
				currentAnnotatedSpliceSiteMap = this.chr2annotatedSpliceSites;
			if(context.isStrandSpecific() && context.getStrand().equals("+"))
				currentAnnotatedSpliceSiteMap = this.chr2annotatedForwardSpliceSites;
			else if(context.isStrandSpecific() && context.getStrand().equals("-"))
				currentAnnotatedSpliceSiteMap = this.chr2annotatedReverseSpliceSites;

			
			//now we check for every full and partial read if it overlaps with an existing split
			int startA = -1;
			int tmpStartA;
			int endA = -1;
			int startB = -1;
			int endB = -1;
			int tmpEndB;
			String readSequence;
			String readSequenceRevComp;
			int mismatches = Integer.MAX_VALUE;
			int minMismatches = -1;
			int insertionSize = -1;
			int minMismatchesOfFullAlignment = -1;
			boolean hasFullReadAlignment = false;
			boolean hasSplicedAlignment = false;
			char hasSpliceSignal = '0';
			boolean foundValidSplit = false;
			boolean hasInsertion = false;
			boolean addedCandidate;
			boolean foundLocationWithKnownSpliceSignal = false;
			boolean foundLocationOverlappingKnownJunction = false;
			
			ArrayList<Integer> splitPositions = new ArrayList<Integer>();
			ArrayList<ReadLocation> fullLocationsToAdd = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> locationsToAdd = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> clippedLocations = new ArrayList<ReadLocation>();
			NavigableMap<Integer,HashSet<Integer>> subMap;
			StringBuilder referenceSequenceBuffer = new StringBuilder();
			StringBuilder readSequenceBuffer = new StringBuilder();
			StringBuilder donorSite = new StringBuilder();
			StringBuilder acceptorSite = new StringBuilder();
			StringBuilder tmpLineBuilder = new StringBuilder();
			int currentLineLength = 0;
			int[] mismatchesArrayForDownstreamSplits = new int[1];
			int[] mismatchesArrayForUpstreamSplits = new int[1];
			int[] mismatchesArrayForVariableSplitPart = new int[1];
			int[] mismatchArrayForVariableInsertionPart = new int[1];
			ArrayList<int[]> poolForUpstreamInsertions = new ArrayList<int[]>();
			ArrayList<int[]> poolForDownstreamInsertions = new ArrayList<int[]>();
			
			
			int startOfVariableSplitPart;
			
			ArrayList<Pair<Integer,Integer>> currentSplits = new ArrayList<Pair<Integer,Integer>>();
			ArrayList<Pair<Integer,Integer>> upstreamSplits = new ArrayList<Pair<Integer,Integer>>();
			ArrayList<Pair<Integer,Integer>> downstreamSplits = new ArrayList<Pair<Integer,Integer>>();
			PairComparator pairComparator = new PairComparator();
			
			int currentSplitPosA;
			int currentSplitPosB;
			int processedReads = 0;
			ReadLocation tmpLocation;
			ReadLocation locationToHash;
			ArrayList<ReadLocation> tmpLocations;
			String tmpLine;
			HashSet<Integer> insertionSizes = new HashSet<Integer>();
			HashSet<String> alreadyCheckedSplits = new HashSet<String>();
			HashSet<String> newlyAddedSplits = new HashSet<String>();
			HashSet<Integer> indicesToSkip = new HashSet<Integer>();
			StringBuilder splitKey = new StringBuilder();
			StringBuilder extendedSplitKey = new StringBuilder();
			String currentInsertionKey;
			String currentLocationKey;
			
			HashMap<Integer,MutableInt> splitSize2frequency = new HashMap<Integer,MutableInt>();
			HashMap<String,ArrayList<ReadLocation>> alreadyDeterminedLocations = new HashMap<String,ArrayList<ReadLocation>>();
			HashMap<String,ArrayList<ReadLocation>> alreadyDeterminedClippings = new HashMap<String,ArrayList<ReadLocation>>();
			HashSet<String> noValidSplitFound = new HashSet<String>();
			int splitSize;
			
			String splitInfo;
			ArrayList<int[]> mismatchArraysForDownstreamInsertions = new ArrayList<int[]>();
			ArrayList<int[]> mismatchArraysForUpstreamInsertions = new ArrayList<int[]>();
			
			BufferedRandomAccessFile braf = null;
			if(alignmentsFilePath != null) {
				braf = new BufferedRandomAccessFile(new File(alignmentsFilePath),"r",1024);
			}
			
			for(InitialRead read : context.getReads()) {
				
				if(read.getLocations() == null && braf != null) {
					parseReadLocations(context,read,contextOffset,localReference.length(),braf);
				}
			
				
				newlyAddedSplits.clear();
				if(alreadyDeterminedLocations.size() > 20000 || alreadyDeterminedClippings.size() > 20000 || noValidSplitFound.size() > 40000) { 
					alreadyDeterminedLocations.clear();
					alreadyDeterminedClippings.clear();
					noValidSplitFound.clear();
				}
				
				readSequence = null;
			
				
				//those splits were already added at the split extension step
				getAlreadyCheckedSplits(read,alreadyCheckedSplits,splitKey);
				
				
				locationsToAdd.clear();
				fullLocationsToAdd.clear();
				clippedLocations.clear();
				addedCandidate = false;
				foundLocationWithKnownSpliceSignal = false;
				foundLocationOverlappingKnownJunction = false;
				readSequenceRevComp = null;
				
				minMismatches = Integer.MAX_VALUE;
				minMismatchesOfFullAlignment = Integer.MAX_VALUE;
				hasFullReadAlignment = false;
				hasSplicedAlignment = false;
				
				for(InitialReadLocation location : read.getLocations()) {
					if(location.getMappingType() == 'F') {
						if(location.getMissmatches() < minMismatches) {
							minMismatches = location.getMissmatches();
							minMismatchesOfFullAlignment = location.getMissmatches();
						}
						
						hasFullReadAlignment = true;
					}
					else if(location.getMappingType() == 'S') {
						hasSplicedAlignment = true;
					}
				}
				
				
				for(InitialReadLocation location : read.getLocations()) {
					if(location.getMappingType() == 'P') {
						
						if(readSequence == null)
							readSequence = this.string2bitset.decompress(this.read2sequence.get(read.getId()));
						
						tmpLineBuilder.setLength(0);
						currentLocationKey = tmpLineBuilder.append(readSequence).append('_').append(location.getStartA()).append('_').append(location.getEndA()).append('_').append(location.getStrand()).toString();						
						if(alreadyDeterminedLocations.containsKey(currentLocationKey)) {
							tmpLocations = alreadyDeterminedLocations.get(currentLocationKey);
							indicesToSkip.clear();
							ReadLocation l;
							for(int i = 0; i < tmpLocations.size(); i++) {
								l = tmpLocations.get(i);
								
								splitKey.setLength(0);
								splitKey.append(l.getStartA()).append('_').append(l.getEndB());
								if(alreadyCheckedSplits.contains(splitKey.toString())) {
									indicesToSkip.add(i);
									continue;
								}
								
								extendedSplitKey.setLength(0);
								extendedSplitKey.append(l.getStartA()).append('_').append(l.getEndA()).append('_').append(l.getStartB()).append('_').append(l.getEndB());								
								
								if(newlyAddedSplits.contains(extendedSplitKey.toString())) {
									indicesToSkip.add(i);
									continue;
								}
								
								if(l.getMismatches() < minMismatches) {
									minMismatches = l.getMismatches();
									locationsToAdd.clear();
									foundLocationWithKnownSpliceSignal = false;
									foundLocationOverlappingKnownJunction = false;
								}
							}
							
							for(int i = 0; i < tmpLocations.size(); i++) {
								if(indicesToSkip.contains(i))
									continue;
								
								l = tmpLocations.get(i);
								if(l.getMismatches() - minMismatches <= this.maxMissmatchDifference) {
									locationsToAdd.add(l);
									if(l.hasSpliceSignal()) foundLocationWithKnownSpliceSignal = true;
									if(l.overlapsKnownJunction()) foundLocationOverlappingKnownJunction = true;
									
									extendedSplitKey.setLength(0);
									extendedSplitKey.append(l.getStartA()).append('_').append(l.getEndA()).append('_').append(l.getStartB()).append('_').append(l.getEndB());
									newlyAddedSplits.add(extendedSplitKey.toString());
								}
							}
							continue;
						}
						
						
						if(!hasFullReadAlignment && !hasSplicedAlignment && alreadyDeterminedClippings.containsKey(currentLocationKey)) {
							tmpLocations = alreadyDeterminedClippings.get(currentLocationKey);
							for(ReadLocation l : tmpLocations) {
								clippedLocations.add(l);
							}
							continue;
						}

						
						startA = location.getStartA();
						endA = location.getEndA();
						subMap = hashedSplits.subMap(startA,true,endA,true);
						
						currentSplits.clear();
						upstreamSplits.clear();
						downstreamSplits.clear();
						startOfVariableSplitPart = Integer.MIN_VALUE;
						insertionSizes.clear();
						
						for(int splitPositionA : subMap.keySet()) {
							for(int splitPositionB : subMap.get(splitPositionA)) {

								if(noValidSplitFound.contains(currentLocationKey + "_" + splitPositionA + "_" + splitPositionB)) {
									continue;
								}
								
								tmpLineBuilder.setLength(0);
								currentInsertionKey = tmpLineBuilder.append(splitPositionA).append('_').append(splitPositionB).toString();
								if(splitPositionB >= splitPositionA || insertionKeys.contains(currentInsertionKey)) {
									if((splitPositionA - startA + 1) == readSequence.length())
										continue;
									
									if(insertionKeys.contains(currentInsertionKey)) {
										if(splitPositionA > startA + readSequence.length()/2) {
											if(readSequence.length() - ((splitPositionA - startA + 1)) - (splitPositionA - splitPositionB + 1) < 0) 
												continue;
										}
										
										else {
											if(splitPositionA - (startA + splitPositionA - splitPositionB + 1) + 1 < 0) {
												continue;
											}
										}
									}
								}
								
								
								if(splitPositionA - splitPositionB + 1 <= this.maxDelSize && insertionKeys.contains(currentInsertionKey)) {
									insertionSizes.add(splitPositionA - splitPositionB + 1);
								}
								
								if(splitPositionB >= splitPositionA || (splitPositionA - splitPositionB + 1 <= this.maxDelSize && insertionKeys.contains(currentInsertionKey))) {
									downstreamSplits.add(new Pair<Integer,Integer>(splitPositionA,splitPositionB));
								}
								
								else if(!insertionKeys.contains(currentInsertionKey)) {
									upstreamSplits.add(new Pair<Integer,Integer>(splitPositionA,splitPositionB));
								}
							}
						}
						
						Collections.sort(downstreamSplits,pairComparator);
						Collections.sort(upstreamSplits,pairComparator);
						
						currentSplits.addAll(downstreamSplits);
						currentSplits.addAll(upstreamSplits);
											
						if(currentSplits.size() > 0) {
							if(readSequence.length() != mismatchesArrayForDownstreamSplits.length) {
								mismatchesArrayForDownstreamSplits = new int[readSequence.length()];
								mismatchesArrayForUpstreamSplits = new int[readSequence.length()];
								mismatchesArrayForVariableSplitPart = new int[readSequence.length()];
								mismatchArrayForVariableInsertionPart = new int[readSequence.length()];
								poolForUpstreamInsertions.clear();
								poolForDownstreamInsertions.clear();
							}
							
							
							if(readSequenceRevComp == null && location.getStrand() == '-') {
								readSequenceBuffer.setLength(0);
								readSequenceBuffer.append(readSequence);
								readSequenceBuffer.reverse();
								for(int i = 0; i < readSequenceBuffer.length(); i++) {
									readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
								}
								readSequenceRevComp = readSequenceBuffer.toString();
							}
														
							if(!downstreamSplits.isEmpty())
								getMismatchesForSplitPartA(mismatchesArrayForDownstreamSplits,(location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,startA,endA,'0');
							
							if(!upstreamSplits.isEmpty() || !insertionSizes.isEmpty())
								getMismatchesForSplitPartB(mismatchesArrayForUpstreamSplits,(location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,startA,endA,'0');
							
							mismatchArraysForDownstreamInsertions.clear();
							mismatchArraysForUpstreamInsertions.clear();
							if(!insertionSizes.isEmpty()) {
								for(int i = 1; i <= this.maxDelSize; i++) {
									if(insertionSizes.contains(i)) {
										int[] tmpDownstream;
										if(mismatchArraysForDownstreamInsertions.size() < poolForDownstreamInsertions.size()) {
											tmpDownstream = poolForDownstreamInsertions.get(mismatchArraysForDownstreamInsertions.size());
										}
										else {
											tmpDownstream = new int[readSequence.length()];
											poolForDownstreamInsertions.add(tmpDownstream);
										}
										getMismatchesForSplitPartA(tmpDownstream,(location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,startA + i,endA + i,'0');
										mismatchArraysForDownstreamInsertions.add(tmpDownstream);
										
										
										int[] tmpUpstream;
										if(mismatchArraysForUpstreamInsertions.size() < poolForUpstreamInsertions.size()) {
											tmpUpstream = poolForUpstreamInsertions.get(mismatchArraysForUpstreamInsertions.size());
										}
										else {
											tmpUpstream = new int[readSequence.length()];
											poolForUpstreamInsertions.add(tmpUpstream);
										}
										
										getMismatchesForSplitPartB(tmpUpstream,(location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,startA - i,endA - i,'0');
										mismatchArraysForUpstreamInsertions.add(tmpUpstream);
									}
									
									else {
										mismatchArraysForDownstreamInsertions.add(null);
										mismatchArraysForUpstreamInsertions.add(null);
									}
								}
							}
						}
						
						for(int i = 0; i < currentSplits.size(); i++) {
							currentSplitPosA = currentSplits.get(i).getFirst();
							currentSplitPosB = currentSplits.get(i).getSecond();

							
							tmpLineBuilder.setLength(0);
							currentInsertionKey = tmpLineBuilder.append(currentSplitPosA).append('_').append(currentSplitPosB).toString();
							//s2 is downstream of s1
							foundValidSplit = false;
							hasInsertion = false;
							if(currentSplitPosB >= currentSplitPosA || insertionKeys.contains(currentInsertionKey)) {
								startA = location.getStartA();
								endA = currentSplitPosA;
								startB = currentSplitPosB;
								endB = startB + (readSequence.length() - (endA - startA + 1)) - 1;
																
								//check if the intron size != 0
								if(startB - endA == 1 && !insertionKeys.contains(currentInsertionKey)) {
									continue;
								}
								
								
								foundValidSplit = true;								
								hasInsertion = false;
								if(insertionKeys.contains(currentInsertionKey)) {
									
									insertionSize = currentSplitPosA - currentSplitPosB + 1;
									hasInsertion = true;
									
									//if the candidate insertion is located in the first half of the read, we have to correct the calculated start and end positions. 
									if(currentSplitPosA <= startA + readSequence.length()/2) {
										startA += insertionSize;
										endB = startB + (readSequence.length() - (endA - startA + 1)) - 1;
										
										//too large insertion sizes
										if(insertionSize > this.maxDelSize || readSequence.length() - ((endA - startA + 1)) - insertionSize < 0) {
											continue;
										}
										
										mismatches = 0;
										if(endA - startA + 1 > 0)
											mismatches = mismatchArraysForDownstreamInsertions.get(insertionSize - 1)[endA - startA];
										
										if(readSequence.length() - ((endA - startA + 1)) - insertionSize > 0) {
											mismatches += mismatchesArrayForUpstreamSplits[readSequence.length() - ((endA - startA + 1)) - 1 - insertionSize];
										}
									}
									
									else {
										
										//we do not consider insertions at the start or the end
										if(endA - startA + 1 < 0) {
											continue;
										}
																				
										mismatches = 0;
										if(endA - startA + 1 > 0)
											mismatches = mismatchesArrayForDownstreamSplits[endA - startA];
										
										if(readSequence.length() - ((endA - startA + 1)) - insertionSize > 0) {
											mismatches += mismatchArraysForUpstreamInsertions.get(insertionSize - 1)[readSequence.length() - ((endA - startA + 1)) - 1 - insertionSize];
										}
									}
									
								}
								
								
								else {
									if(endB != startOfVariableSplitPart) {
										getMismatchesForSplitPartB(mismatchesArrayForVariableSplitPart,(location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,endB - readSequence.length() + 1,endB,'0');
										startOfVariableSplitPart = endB;
									}
									
									mismatches = mismatchesArrayForDownstreamSplits[endA - startA] + mismatchesArrayForVariableSplitPart[readSequence.length() - ((endA - startA + 1)) - 1];
								}
								
							}
							//s2 is upstream of s1
							else if(currentSplitPosA > location.getStartA()) {
								
								endB = location.getStartA() + readSequence.length() - 1;
								startB = currentSplitPosA;
								endA = currentSplitPosB;
								startA = endA - (readSequence.length() - (endB - startB + 1)) + 1;
								
								//check if the second part consists of the whole read
								if(endB - startB + 1 == readSequence.length())
									continue;
								
								
								//check if the intron size != 0
								if(startB - endA == 1)
									continue;
								
								
								foundValidSplit = true;
								if(startA != startOfVariableSplitPart) {
									getMismatchesForSplitPartA(mismatchesArrayForVariableSplitPart,(location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,startA,startA + readSequence.length() - 1,'0');
									startOfVariableSplitPart = startA;
								}
								
								mismatches = mismatchesArrayForVariableSplitPart[endA - startA] + mismatchesArrayForUpstreamSplits[readSequence.length() - ((endA - startA + 1)) - 1];
							}
							
							
							if(foundValidSplit) {
								splitKey.setLength(0);
								splitKey.append(startA);
								splitKey.append('_');
								splitKey.append(endB);

								extendedSplitKey.setLength(0);
								extendedSplitKey.append(startA).append('_').append(endA).append('_').append(startB).append('_').append(endB);								
							}
							
							if(foundValidSplit  && !alreadyCheckedSplits.contains(splitKey.toString()) && !newlyAddedSplits.contains(extendedSplitKey.toString())) {
								//newlyAddedSplits.add(new int[]{startA,endA,startB,endB});
								newlyAddedSplits.add(extendedSplitKey.toString());
								if(mismatches <= this.maxMissmatches) {
									
									//valid hit with too many missmatches
									if(mismatches > minMismatches && (mismatches - minMismatches) > this.maxMissmatchDifference) {
										
										//hit could be an allowed hit for the next read
										if(!alreadyDeterminedLocations.containsKey(currentLocationKey))
											alreadyDeterminedLocations.put(currentLocationKey,new ArrayList<ReadLocation>());
										
										locationToHash = new ReadLocation(context.getChr(), location.getStrand(),'S',startA,endA,startB,endB,mismatches);
										if(!hasInsertion && (startB - endA - 1) > this.maxDelSize) {
											hasSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA,startB,locationToHash.getStrand(),context.isStrandSpecific());
											locationToHash.setStrandOfSpliceSignal(hasSpliceSignal);
											
											if(currentAnnotatedSpliceSiteMap != null && currentAnnotatedSpliceSiteMap.containsKey(locationToHash.getChr())) {
												if((currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).containsKey(locationToHash.getEndA()) && currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).get(locationToHash.getEndA()).contains(locationToHash.getStartB())) ||
												   (currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).containsKey(locationToHash.getStartB()) && currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).get(locationToHash.getStartB()).contains(locationToHash.getEndA()))){
													locationToHash.setOverlapsKnownJunction(true);
												}
											}
										}
										locationToHash.setOverlapsKnownJunction(locationToHash.overlapsKnownJunction());
										locationToHash.setStrandOfSpliceSignal(locationToHash.getStrandOfSpliceSignal());
										alreadyDeterminedLocations.get(currentLocationKey).add(locationToHash);
										
										continue;
									}
									
									//valid hit with new min missmatch boundary
									if(mismatches < minMismatches && (minMismatches - mismatches) > this.maxMissmatchDifference) {
										minMismatches = mismatches;
										locationsToAdd.clear();
										fullLocationsToAdd.clear();
										foundLocationWithKnownSpliceSignal = false;
										foundLocationOverlappingKnownJunction = false;
									}
									
									if(locationsToAdd.size() < readLocationPool.size()) {
										tmpLocation = readLocationPool.get(locationsToAdd.size());
										tmpLocation.updateLocation(context.getChr(), location.getStrand(),'S',startA,endA,startB,endB,mismatches);
									}
									else {
										tmpLocation = new ReadLocation(context.getChr(), location.getStrand(),'S',startA,endA,startB,endB,mismatches);
										readLocationPool.add(tmpLocation);
									}
									
									if(!hasInsertion && (startB - endA - 1) > this.maxDelSize) {
										hasSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA,startB,tmpLocation.getStrand(),context.isStrandSpecific());
										tmpLocation.setStrandOfSpliceSignal(hasSpliceSignal);
										if(tmpLocation.hasSpliceSignal())
											foundLocationWithKnownSpliceSignal = true;
										
										//in case an annotation is given check if the actual split overlaps with an annotated junction
										if(currentAnnotatedSpliceSiteMap != null && currentAnnotatedSpliceSiteMap.containsKey(tmpLocation.getChr())) {
											if((currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getEndA()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getEndA()).contains(tmpLocation.getStartB())) ||
											   (currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getStartB()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getStartB()).contains(tmpLocation.getEndA()))){
												tmpLocation.setOverlapsKnownJunction(true);
												foundLocationOverlappingKnownJunction = true;
											}
										}
									}
									
									locationsToAdd.add(tmpLocation);
									if(!alreadyDeterminedLocations.containsKey(currentLocationKey))
										alreadyDeterminedLocations.put(currentLocationKey,new ArrayList<ReadLocation>());
									
									locationToHash = new ReadLocation(context.getChr(), tmpLocation.getStrand(),tmpLocation.getMappingType(),tmpLocation.getStartA(),tmpLocation.getEndA(),tmpLocation.getStartB(),tmpLocation.getEndB(),tmpLocation.getMismatches());
									locationToHash.setOverlapsKnownJunction(tmpLocation.overlapsKnownJunction());
									locationToHash.setStrandOfSpliceSignal(tmpLocation.getStrandOfSpliceSignal());
									alreadyDeterminedLocations.get(currentLocationKey).add(locationToHash);
								}
							}
							
							else if(foundValidSplit && !newlyAddedSplits.contains(extendedSplitKey.toString())) {
								if(mismatches <= this.maxMissmatches) {
									//hit could be an allowed hit for the next read
									if(!alreadyDeterminedLocations.containsKey(currentLocationKey))
										alreadyDeterminedLocations.put(currentLocationKey,new ArrayList<ReadLocation>());
									
									locationToHash = new ReadLocation(context.getChr(), location.getStrand(),'S',startA,endA,startB,endB,mismatches);
									if(!hasInsertion && (startB - endA - 1) > this.maxDelSize) {
										hasSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA,startB,locationToHash.getStrand(),context.isStrandSpecific());
										locationToHash.setStrandOfSpliceSignal(hasSpliceSignal);
										
										if(currentAnnotatedSpliceSiteMap != null && currentAnnotatedSpliceSiteMap.containsKey(locationToHash.getChr())) {
											if((currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).containsKey(locationToHash.getEndA()) && currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).get(locationToHash.getEndA()).contains(locationToHash.getStartB())) ||
											   (currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).containsKey(locationToHash.getStartB()) && currentAnnotatedSpliceSiteMap.get(locationToHash.getChr()).get(locationToHash.getStartB()).contains(locationToHash.getEndA()))){
												locationToHash.setOverlapsKnownJunction(true);
											}
										}
									}
									locationToHash.setOverlapsKnownJunction(locationToHash.overlapsKnownJunction());
									locationToHash.setStrandOfSpliceSignal(locationToHash.getStrandOfSpliceSignal());
									alreadyDeterminedLocations.get(currentLocationKey).add(locationToHash);
								}
							}
							
							
							else if(!foundValidSplit) {
								noValidSplitFound.add(currentLocationKey + "_" + currentSplitPosA + "_" + currentSplitPosB);
							}
							
						}
						
						
						if(this.clipping && locationsToAdd.isEmpty() && !hasFullReadAlignment && !hasSplicedAlignment && !read.getId().contains("MSC")) {
							
							if(readSequenceRevComp == null && location.getStrand() == '-') {
								readSequenceBuffer.setLength(0);
								readSequenceBuffer.append(readSequence);
								readSequenceBuffer.reverse();
								for(int i = 0; i < readSequenceBuffer.length(); i++) {
									readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
								}
								readSequenceRevComp = readSequenceBuffer.toString();
							}
							
							ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices =  getBestLocalAlignmentIndices((location.getStrand() == '+')? readSequence : readSequenceRevComp,localReference,contextOffset, readSequenceBuffer, location.getStartA(), location.getStartA() + readSequence.length() -1 ,'0', 1, -4, this.maxMissmatches,this.seedLength,false,false);
							for(int j = 0; j < clippingIndices.size(); j++) {
								clippedLocations.add(new ReadLocation(context.getChr(), location.getStrand(),'F',location.getStartA(),location.getStartA() + readSequence.length() -1,-1 * clippingIndices.get(j).getFirst(),-1 * clippingIndices.get(j).getSecond(),clippingIndices.get(j).getThird()));
								
								
								if(!alreadyDeterminedClippings.containsKey(currentLocationKey))
									alreadyDeterminedClippings.put(currentLocationKey,new ArrayList<ReadLocation>());
								alreadyDeterminedClippings.get(currentLocationKey).add(new ReadLocation(context.getChr(), location.getStrand(),'F',location.getStartA(),location.getStartA() + readSequence.length() -1,-1 * clippingIndices.get(j).getFirst(),-1 * clippingIndices.get(j).getSecond(),clippingIndices.get(j).getThird()));
							}

						}
					}
				}
				
				
				if((hasFullReadAlignment && locationsToAdd.size() >= 2 && minMismatches >= minMismatchesOfFullAlignment) ||
						   (!hasFullReadAlignment && locationsToAdd.size() > 2)) {
							splitSize2frequency.clear();
							for(ReadLocation l : locationsToAdd) {
								splitSize = l.getEndA() - l.getStartA() + 1;
								if(splitSize2frequency.containsKey(splitSize))
									splitSize2frequency.get(splitSize).increment();
								else
									splitSize2frequency.put(splitSize, new MutableInt(1));
							}
							
							for(MutableInt freq : splitSize2frequency.values()) {
								if(freq.intValue() > locationsToAdd.size()/2) {
									locationsToAdd.clear();
									break;
								}
							}
						}
				
				if(!locationsToAdd.isEmpty()) {
					context.getPartialAndFullReads2filePointer().put(read.getId(), filePointer);
				}
				
				else if(this.clipping && !clippedLocations.isEmpty() && locationsToAdd.isEmpty()) {
					context.getPartialAndFullReads2filePointer().put(read.getId(), filePointer);
					
					tmpLineBuilder.setLength(0);
					tmpLineBuilder.append(context.getId()).append("\t").append(read.getId()).append("\t").append(context.getChr()).append("\t");
					currentLineLength = tmpLineBuilder.length();
					for(ReadLocation locationToAdd : clippedLocations) {
						tmpLineBuilder.setLength(currentLineLength);
						tmpLine = tmpLineBuilder.append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
						pw.write(tmpLine);
						pw.newLine();
						filePointer += tmpLine.length() + 1;
					}
				}
				
				
				//print out new locations
				if(!locationsToAdd.isEmpty()) {
					tmpLineBuilder.setLength(0);
					tmpLineBuilder.append(context.getId()).append("\t").append(read.getId()).append("\t").append(context.getChr()).append("\t");
					currentLineLength = tmpLineBuilder.length();
				}
				
				for(ReadLocation locationToAdd : locationsToAdd) {
					if(currentAnnotatedSpliceSiteMap == null || !foundLocationOverlappingKnownJunction) {
						if(!this.preferExtensionsWithKnownSpliceSignal || !foundLocationWithKnownSpliceSignal || locationToAdd.hasSpliceSignal()) {
							tmpLineBuilder.setLength(currentLineLength);
							tmpLine = tmpLineBuilder.append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
							pw.write(tmpLine);
							pw.newLine();
							filePointer += tmpLine.length() + 1;
						}
					}
					
					else if(locationToAdd.overlapsKnownJunction()) {
						tmpLineBuilder.setLength(currentLineLength);
						tmpLine = tmpLineBuilder.append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
						pw.write(tmpLine);
						pw.newLine();
						filePointer += tmpLine.length() + 1;
					}
				}
				
				
				
				processedReads++;
				
				
				//print out new locations for duplicated reads
				if(read.getDuplicates() != null) {
					for(String id : read.getDuplicates()) {
						if(!locationsToAdd.isEmpty()) {
							context.getPartialAndFullReads2filePointer().put(id, filePointer);
						}
						
						else if(this.clipping && !clippedLocations.isEmpty() && locationsToAdd.isEmpty()) {
							context.getPartialAndFullReads2filePointer().put(id, filePointer);
							
							tmpLineBuilder.setLength(0);
							tmpLineBuilder.append(context.getId()).append("\t").append(id).append("\t").append(context.getChr()).append("\t");
							currentLineLength = tmpLineBuilder.length();
							for(ReadLocation locationToAdd : clippedLocations) {
								tmpLineBuilder.setLength(currentLineLength);
								tmpLine = tmpLineBuilder.append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
								pw.write(tmpLine);
								pw.newLine();
								filePointer += tmpLine.length() + 1;
							}
						}
						
						
						if(!locationsToAdd.isEmpty()) {
							tmpLineBuilder.setLength(0);
							tmpLineBuilder.append(context.getId()).append("\t").append(id).append("\t").append(context.getChr()).append("\t");
							currentLineLength = tmpLineBuilder.length();
						}
						
						for(ReadLocation locationToAdd : locationsToAdd) {
							if(currentAnnotatedSpliceSiteMap == null || !foundLocationOverlappingKnownJunction) {
								if(!this.preferExtensionsWithKnownSpliceSignal || !foundLocationWithKnownSpliceSignal || locationToAdd.hasSpliceSignal()) {
									tmpLineBuilder.setLength(currentLineLength);
									tmpLine = tmpLineBuilder.append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString(); 
									pw.write(tmpLine);
									pw.newLine();
									filePointer += tmpLine.length() + 1;
								}
							}
							
							else if(locationToAdd.overlapsKnownJunction()) {
								tmpLineBuilder.setLength(0);
								tmpLine = tmpLineBuilder.append(context.getId()).append("\t").append(id).append("\t").append(locationToAdd.getChr()).append("\t").append(locationToAdd.getStrand()).append("\t").append(locationToAdd.getStartA()).append("\t").append(locationToAdd.getEndA()).append("\t").append(locationToAdd.getStartB()).append("\t").append(locationToAdd.getEndB()).append("\t").append(locationToAdd.getMismatches()).append("\t").append(locationToAdd.getStrandOfSpliceSignal()).append("\t").append(locationToAdd.overlapsKnownJunction()).append("\t").append(read.getOverallMappingCount()).append("\t").append(read.getOverallValidPairsCount()).toString();
								pw.write(tmpLine);
								pw.newLine();
								filePointer += tmpLine.length() + 1;
							}
						}
					}
				}
				
				//free mem
				for(InitialReadLocation location : read.getLocations()) {
					location = null;
				}
				read.getLocations().clear();
				if(read.getDuplicates() != null) 
					read.getDuplicates().clear();
				
			}
			
			if(braf != null) braf.close();
			
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	
	private class PairComparator implements Comparator {

		@Override
		public int compare(Object o1, Object o2) {
			Pair<Integer,Integer> p1 = (Pair<Integer,Integer>)o1;
			Pair<Integer,Integer> p2 = (Pair<Integer,Integer>)o2;
			return p1.getSecond().compareTo(p2.getSecond());
		}
		
	}
	

	
	private ArrayList<ReadLocation> getSplitLocations(int start, int end, String readSequence, StringBuilder localReference, int[] mismatchesArrayForDownstreamSplits, int[] mismatchesArrayForUpstreamSplits, int contextOffset, String chr, char strand, TreeMap<Integer,HashSet<Integer>> hashedSplits, HashSet<String> insertionKeys, HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentAnnotatedSpliceSiteMap,ArrayList<ReadLocation> readLocationPool, ArrayList<int[]> mismatchArrayPool, boolean strandSpecific,boolean onlyUpstream, boolean onlyDownstream) {
		//now we check for every full and partial read if it overlaps with an existing split
		
		String readSequenceRevComp = null;
		
		int startA = -1;
		int tmpStartA;
		int endA = -1;
		int startB = -1;
		int endB = -1;		
		int tmpEndB;
		
		int insertionSize;
		int mismatches = Integer.MAX_VALUE;
		int minMismatches = Integer.MAX_VALUE;
		boolean hasInsertion = false;
		char hasSpliceSignal = '0';
		boolean foundValidSplit = false;
		
		ArrayList<Integer> splitPositions = new ArrayList<Integer>();
		ArrayList<ReadLocation> locationsToAdd = new ArrayList<ReadLocation>();
		NavigableMap<Integer,HashSet<Integer>> subMap;
		
		StringBuilder readSequenceBuffer = new StringBuilder();
		StringBuilder donorSite = new StringBuilder();
		StringBuilder acceptorSite = new StringBuilder();
		if(mismatchesArrayForDownstreamSplits.length != readSequence.length())
			mismatchesArrayForDownstreamSplits = new int[readSequence.length()];
		
		if(mismatchesArrayForUpstreamSplits.length != readSequence.length())
			mismatchesArrayForUpstreamSplits = new int[readSequence.length()];
		
		int[] mismatchesArrayForVariableSplitPart = new int[readSequence.length()];
		
		ArrayList<int[]> mismatchArraysForDownstreamInsertions = new ArrayList<int[]>();
		if(!mismatchArrayPool.isEmpty() && mismatchArrayPool.get(0).length != readSequence.length())
			mismatchArrayPool.clear();
		
		
		int startOfVariableSplitPart;
		ArrayList<Split> currentSplits = new ArrayList<Split>();
		int currentSplitPosA;
		int currentSplitPosB;
		ReadLocation tmpLocation;
		HashSet<String> alreadyCheckedSplits = new HashSet<String>();
		HashSet<Integer> insertionSizes = new HashSet<Integer>();
		HashSet<String> newlyAddedSplits = new HashSet<String>();
		StringBuilder splitKey = new StringBuilder();
		String currentInsertionKey;
		StringBuilder extendedSplitKey = new StringBuilder();
		
		
		subMap = hashedSplits.subMap(start,true,end,true);
		currentSplits.clear();
		startOfVariableSplitPart = Integer.MIN_VALUE;
		for(int splitPositionA : subMap.keySet()) {
			for(int splitPositionB : subMap.get(splitPositionA)) {
				
				splitKey.setLength(0);
				splitKey.append(splitPositionA).append('_').append(splitPositionB);
				currentInsertionKey = splitKey.toString();
				if((!onlyUpstream && !onlyDownstream)  || 
					(onlyUpstream &&  (splitPositionB < splitPositionA)) ||
					(onlyDownstream && splitPositionB >= splitPositionA) ||
					insertionKeys.contains(currentInsertionKey)) {
					
					if(splitPositionB >= splitPositionA || insertionKeys.contains(currentInsertionKey)) {
						if((splitPositionA - start + 1) == readSequence.length())
							continue;
						
						if(insertionKeys.contains(currentInsertionKey)) {
							if(splitPositionA > start + readSequence.length()/2) {
								if(readSequence.length() - ((splitPositionA - start + 1)) - (splitPositionA - splitPositionB + 1) < 0) 
									continue;
							}
							
							else {
								if(splitPositionA - (start + splitPositionA - splitPositionB + 1) + 1 < 0) {
									continue;
								}
							}
						}
					}
					
					currentSplits.add(new Split(splitPositionA,splitPositionB));
					if(insertionKeys.contains(currentInsertionKey)) {
						insertionSizes.add(splitPositionA - splitPositionB + 1);
					}
				}
			}
		}
							
		if(currentSplits.size() > 0) {
			Collections.sort(currentSplits);
			if(strand == '-') {
				readSequenceBuffer.setLength(0);
				readSequenceBuffer.append(readSequence);
				readSequenceBuffer.reverse();
				for(int i = 0; i < readSequenceBuffer.length(); i++) {
					readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
				}
				readSequenceRevComp = readSequenceBuffer.toString();
			}
			
			getMismatchesForSplitPartA(mismatchesArrayForDownstreamSplits,(strand == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,start,end,'0');
			getMismatchesForSplitPartB(mismatchesArrayForUpstreamSplits,(strand == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,start,end,'0');
			
			mismatchArraysForDownstreamInsertions.clear();
			for(int i = 1; i <= this.maxDelSize; i++) {
				if(insertionSizes.contains(i)) {
					int[] tmpDownstream;
					if(mismatchArraysForDownstreamInsertions.size() < mismatchArrayPool.size())
						tmpDownstream = mismatchArrayPool.get(mismatchArraysForDownstreamInsertions.size());
					
					else {
						tmpDownstream = new int[readSequence.length()];
						mismatchArrayPool.add(tmpDownstream);
						}
					
					getMismatchesForSplitPartA(tmpDownstream,(strand == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,start + i,end + i,'0');
					mismatchArraysForDownstreamInsertions.add(tmpDownstream);
				}
				else {
					mismatchArraysForDownstreamInsertions.add(null);
				}
			}
			
		}
		
		for(int i = 0; i < currentSplits.size(); i++) {
			currentSplitPosA = currentSplits.get(i).getSplitPositionA();
			currentSplitPosB = currentSplits.get(i).getSplitPositionB();
			
			splitKey.setLength(0);
			splitKey.append(currentSplitPosA).append('_').append(currentSplitPosB);
			currentInsertionKey = splitKey.toString();
			
				//s2 is downstream of s1
				foundValidSplit = false;
				hasInsertion = false;
				if(currentSplitPosB >= currentSplitPosA || insertionKeys.contains(currentInsertionKey)) {
					startA = start;
					endA = currentSplitPosA;
					startB = currentSplitPosB;
					endB = startB + (readSequence.length() - (endA - startA + 1)) - 1;

					foundValidSplit = true;
					
					if(insertionKeys.contains(currentInsertionKey)) {										
						insertionSize = currentSplitPosA - currentSplitPosB + 1;
						hasInsertion = true;
					
					    //if the candidate insertion is located in the first half of the read, we have to correct the calculated start and end positions. 
					    if(currentSplitPosA <= startA + readSequence.length()/2) {
					    	startA += insertionSize;
					    	endB = startB + (readSequence.length() - (endA - startA + 1)) - 1;
						
						    //too large insertion sizes
						    if(insertionSize > this.maxDelSize || readSequence.length() - ((endA - startA + 1)) - insertionSize < 0) {
								continue;
							}
							
						    mismatches = 0;
						    if(endA - startA + 1 > 0)
						    	mismatches = mismatchArraysForDownstreamInsertions.get(insertionSize - 1)[endA - startA];
							if(readSequence.length() - ((endA - startA + 1)) - insertionSize > 0) {
								mismatches += mismatchesArrayForUpstreamSplits[readSequence.length() - ((endA - startA + 1)) - 1 - insertionSize];
							}
					    }
					
					    else {
						
							//we do not consider insertions at the start or the end
					    	if(endA - startA + 1 < 0 || insertionSize > this.maxDelSize)
								continue;
					
							mismatches = 0;
							if(endA - startA + 1 > 0)
								mismatches = mismatchesArrayForDownstreamSplits[endA - startA];
							if(readSequence.length() - ((endA - startA + 1)) - insertionSize > 0) {
								if(endB != startOfVariableSplitPart) {
									getMismatchesForSplitPartB(mismatchesArrayForVariableSplitPart,(strand == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,endB - readSequence.length() + 1,endB,'0');
									startOfVariableSplitPart = endB;
								}
								mismatches += mismatchesArrayForVariableSplitPart[readSequence.length() - ((endA - startA + 1)) - 1 - insertionSize];
							}
					    }
						
				    }
				
					else {
						if(endB != startOfVariableSplitPart) {
							getMismatchesForSplitPartB(mismatchesArrayForVariableSplitPart,(strand == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,endB - readSequence.length() + 1,endB,'0');
							startOfVariableSplitPart = endB;
						}
						mismatches = mismatchesArrayForDownstreamSplits[endA - startA] + mismatchesArrayForVariableSplitPart[readSequence.length() - ((endA - startA + 1)) - 1];
					}
				}
				
				//s2 is upstream of s1
				else if(currentSplitPosA > start) {
					hasInsertion = false;
					endB = start + readSequence.length() - 1;
					startB = currentSplitPosA;
					endA = currentSplitPosB;
					startA = endA - (readSequence.length() - (endB - startB + 1)) + 1;
					
					if(endB - startB + 1 == readSequence.length())
						continue;
					
					
					foundValidSplit = true;
					
					if(startA != startOfVariableSplitPart) {
						getMismatchesForSplitPartA(mismatchesArrayForVariableSplitPart,(strand == '+')? readSequence : readSequenceRevComp,localReference,contextOffset,null,startA,startA + readSequence.length() - 1,'0');
						startOfVariableSplitPart = startA;
					}
					
				
					mismatches = mismatchesArrayForVariableSplitPart[endA - startA] + mismatchesArrayForUpstreamSplits[readSequence.length() - ((endA - startA + 1)) - 1];
				}
				
				
				splitKey.setLength(0);
				splitKey.append(startA);
				splitKey.append('_');
				splitKey.append(endB);
				
				extendedSplitKey.setLength(0);
				extendedSplitKey.append(startA).append('_').append(endA).append('_').append(startB).append('_').append(endB);
				
				
				if(foundValidSplit  && !alreadyCheckedSplits.contains(splitKey.toString()) && !newlyAddedSplits.contains(extendedSplitKey.toString())) {
					newlyAddedSplits.add(extendedSplitKey.toString());
					if(mismatches <= this.maxMissmatches) {
						
						//valid hit with too many missmatches
						if(mismatches > minMismatches && (mismatches - minMismatches) > this.maxMissmatchDifference)
							continue;
						
						//valid hit with new min missmatch boundary
						if(mismatches < minMismatches && (minMismatches - mismatches) > this.maxMissmatchDifference) {
							minMismatches = mismatches;
							locationsToAdd.clear();
						}
						
						if(locationsToAdd.size() < readLocationPool.size()) {
							tmpLocation = readLocationPool.get(locationsToAdd.size());
							tmpLocation.updateLocation(chr, strand,'S',startA,endA,startB,endB,mismatches);
						}
						else {
							tmpLocation = new ReadLocation(chr, strand,'S',startA,endA,startB,endB,mismatches);
							readLocationPool.add(tmpLocation);
						}
						
						if(!hasInsertion && (startB - endA - 1) > this.maxDelSize) {
							hasSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA,startB,tmpLocation.getStrand(),strandSpecific);
							tmpLocation.setStrandOfSpliceSignal(hasSpliceSignal);
							
							
							//in case an annotation is given check if the actual split overlaps with an annotated junction
							if(currentAnnotatedSpliceSiteMap != null && currentAnnotatedSpliceSiteMap.containsKey(tmpLocation.getChr())) {
								if((currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getEndA()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getEndA()).contains(tmpLocation.getStartB())) ||
								   (currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).containsKey(tmpLocation.getStartB()) && currentAnnotatedSpliceSiteMap.get(tmpLocation.getChr()).get(tmpLocation.getStartB()).contains(tmpLocation.getEndA()))){
									tmpLocation.setOverlapsKnownJunction(true);
								}
							}
						}
						
						locationsToAdd.add(tmpLocation);												
					}
				}
		}
	
		return locationsToAdd;
	}

		
	public void extendSplitAndFullAlignments(String multiMappingFilePath, Context context, StringBuilder localReference, int contextOffset, String outputFilePath, int minStartPositionOverlaps) {
		try {
			
			//hash splits
			TreeMap<Integer,HashSet<Integer>> hashedSplits = new TreeMap<Integer,HashSet<Integer>>();
			HashSet<String> insertionKeys = new HashSet<String>();
			HashMap<String,Split> key2split = new HashMap<String,Split>();
			String splitKey;
			String currentLine;
			String[] splittedLine;
			StringBuffer sb = new StringBuffer();
			Pattern tabPattern = Pattern.compile("\t");
			BufferedReader br = new BufferedReader(new FileReader(new File(multiMappingFilePath)));
			
			int endA;
			int startB;
			char strandOfSpliceSignal;
			boolean overlapsKnownJunction;
			
			StringBuilder tmpBuilder = new StringBuilder();
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				endA= Integer.valueOf(splittedLine[5]);
				startB = Integer.valueOf(splittedLine[6]);
				
				//full read
				if(startB <= 0)
					continue;
				
				strandOfSpliceSignal = splittedLine[9].charAt(0);
				overlapsKnownJunction = Boolean.valueOf(splittedLine[10]);
				//no signal, check if it is an indel, otherwise skip this split
				if(strandOfSpliceSignal == '0' && !overlapsKnownJunction) {
					if(startB > endA && (startB - endA) > this.maxDelSize)
						continue;
				}
				
				
				tmpBuilder.setLength(0);
				splitKey = tmpBuilder.append(splittedLine[5]).append('_').append(splittedLine[6]).toString();
				if(!key2split.containsKey(splitKey)) {
					key2split.put(splitKey, new Split(Integer.valueOf(splittedLine[5]),Integer.valueOf(splittedLine[6])));
				}
				key2split.get(splitKey).addStartPosition(Integer.valueOf(splittedLine[4]));
			}
			br.close();
			
			for(Split split : key2split.values()) {
				if(split.getDifferentStartPositions() < minStartPositionOverlaps)
					continue;
					
					
				if(hashedSplits.containsKey(split.getSplitPositionA()))
					hashedSplits.get(split.getSplitPositionA()).add(split.getSplitPositionB());
				else {
					HashSet<Integer> tmpSet = new HashSet<Integer>();
					tmpSet.add(split.getSplitPositionB());
					hashedSplits.put(split.getSplitPositionA(), tmpSet);
				}
				
				
				if(split.getSplitPositionA() < split.getSplitPositionB()) {
					if(hashedSplits.containsKey(split.getSplitPositionB())) {
						hashedSplits.get(split.getSplitPositionB()).add(split.getSplitPositionA());
					}
					else {
						HashSet<Integer> tmpSet = new HashSet<Integer>();
						tmpSet.add(split.getSplitPositionA());
						hashedSplits.put(split.getSplitPositionB(), tmpSet);
					}
				}
				
				else {
					tmpBuilder.setLength(0);
					tmpBuilder.append(split.getSplitPositionA()).append('_').append(split.getSplitPositionB());
					insertionKeys.add(tmpBuilder.toString());
				}
				
			}
			key2split.clear();

			
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentAnnotatedSpliceSiteMap = null;
			if(!context.isStrandSpecific())
				currentAnnotatedSpliceSiteMap = this.chr2annotatedSpliceSites;
			if(context.isStrandSpecific() && context.getStrand().equals("+"))
				currentAnnotatedSpliceSiteMap = this.chr2annotatedForwardSpliceSites;
			else if(context.isStrandSpecific() && context.getStrand().equals("-"))
				currentAnnotatedSpliceSiteMap = this.chr2annotatedReverseSpliceSites;

			String readSequence;
			String currentFullId;
			String prevFullId = null;
			
			String currentSplitId;
			String prevSplitId = null;
			String[] splittedLocation;
			
			boolean foundLocationWithKnownSpliceSignal = false;
			boolean foundLocationOverlappingKnownJunction = false;
			
			ArrayList<ReadLocation> locationsToAdd = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> prevFullLocationsToAdd = new ArrayList<ReadLocation>();
			ArrayList<String> prevSplitLocationsToAdd = new ArrayList<String>();
			ArrayList<ReadLocation> locationsToRemove = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> locationsToAddA = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> locationsToAddB = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> currentCombinations = new ArrayList<ReadLocation>();
			ReadLocation tmpLocation;
			int tmpMismatches = Integer.MAX_VALUE;
			int minMismatches;
			int mismatchesOfOriginalAlignment;
			int tmpMinMismatches;
			int minMismatchesOfSplits;
			
			
			StringBuilder donorSite = new StringBuilder();
			StringBuilder acceptorSite = new StringBuilder();
			
			String chr;
			char strand;
			int startA;
			int endB;
			int shift;
			boolean artificialSplit;
			char tmpStrandOfSpliceSignal;
			
			ArrayList<ReadLocation> readLocationPool = new ArrayList<ReadLocation>();
			ArrayList<int[]> mismatchArrayPool = new ArrayList<int[]>();
			int[] mismatchesArrayForUpstreamSplits = new int[1];
			int[] mismatchesArrayForDownstreamSplits = new int[1];
			
			
			br = new BufferedReader(new FileReader(new File(multiMappingFilePath)));
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outputFilePath),true)), 1024 * 1024 * 10);
			HashMap<Integer,MutableInt> splitSize2frequency = new HashMap<Integer,MutableInt>();
			int splitSize;
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				
				//multi split candidate (processed later)
				if(splittedLine[1].contains("::MSC::"))
					continue;
			
				
				chr = splittedLine[2];
				strand = splittedLine[3].charAt(0);
				startA = Integer.valueOf(splittedLine[4]);
				endA = Integer.valueOf(splittedLine[5]);
				startB = Integer.valueOf(splittedLine[6]);
				endB = Integer.valueOf(splittedLine[7]);
				
				//clipped read
				if(startB < 0 || endB < 0)
					continue;
				
				minMismatches = Integer.valueOf(splittedLine[8]);
				mismatchesOfOriginalAlignment = Integer.valueOf(splittedLine[8]);
				strandOfSpliceSignal = splittedLine[9].charAt(0);
				overlapsKnownJunction = Boolean.valueOf(splittedLine[10]);
				readSequence = this.string2bitset.decompress(this.read2sequence.get(splittedLine[1]));
				
				tmpBuilder.setLength(0);
				tmpBuilder.append(readSequence).append('_').append(startA).append('_').append(endA);
				currentFullId = tmpBuilder.toString();
				
				tmpBuilder.append('_').append(startB).append('_').append(endB);
				currentSplitId = tmpBuilder.toString();
				
				//processing full read
				if(startB == 0) {
					if(prevFullId != null && currentFullId.equals(prevFullId)) {
						for(ReadLocation location : prevFullLocationsToAdd) {
							sb.setLength(0);
							sb.append(context.getId()).append("\t").append(splittedLine[1]).append("\t").append(chr).append("\t").append(strand).append("\t").append(location.getStartA()).append(",").append(location.getEndA()).append(",").append(location.getStartB()).append(",").append(location.getEndB()).append("\t").append(location.getMismatches()).append("\t").append(location.getStrandOfSpliceSignal()).append("\t").append(location.overlapsKnownJunction()).append("\t").append(splittedLine[11]).append("\t").append(splittedLine[12]).append("\t").append(splittedLine[13]).append("\n");
							bw.write(sb.toString());
						}
						continue;
					}
					
					locationsToAddA = getSplitLocations(startA, endA, readSequence,localReference,mismatchesArrayForDownstreamSplits,mismatchesArrayForUpstreamSplits,contextOffset,chr, strand, hashedSplits,insertionKeys, currentAnnotatedSpliceSiteMap,readLocationPool,mismatchArrayPool, context.isStrandSpecific(),false,false);
					
					locationsToAdd.clear();
					locationsToRemove.clear();
					foundLocationWithKnownSpliceSignal = false;
					foundLocationOverlappingKnownJunction = false;
					minMismatchesOfSplits = Integer.MAX_VALUE;
					for(ReadLocation location : locationsToAddA) {
						if(location.getMismatches() > mismatchesOfOriginalAlignment + this.maxMissmatchDifference) {
							locationsToRemove.add(location);
							continue;
						}
						
						if(location.getMismatches() < minMismatchesOfSplits)
							minMismatchesOfSplits = location.getMismatches();
						
						if(location.hasSpliceSignal())
							foundLocationWithKnownSpliceSignal = true;
						
						if(location.overlapsKnownJunction())
							foundLocationOverlappingKnownJunction = true;
						
					}
					locationsToAddA.removeAll(locationsToRemove);
					
					
					if(locationsToAddA.size() >= 2) {
						if(minMismatchesOfSplits >= mismatchesOfOriginalAlignment) {
							splitSize2frequency.clear();
							for(ReadLocation l : locationsToAddA) {
								splitSize = l.getEndA() - l.getStartA() + 1;
								if(splitSize2frequency.containsKey(splitSize))
									splitSize2frequency.get(splitSize).increment();
								else
									splitSize2frequency.put(splitSize, new MutableInt(1));
							}
							
							for(MutableInt freq : splitSize2frequency.values()) {
								if(freq.intValue() > locationsToAddA.size()/2) {
									locationsToAddA.clear();
									break;
								}
							}
						}
					}
					
					
					//first: if the original alignment has 0 mismatches, we only consider candidates overlapping known splice signals or annotated splice sites
					if(mismatchesOfOriginalAlignment == 0) {
						locationsToRemove.clear();
						for(ReadLocation location : locationsToAddA) {
							if(!location.overlapsKnownJunction() && !location.hasSpliceSignal())
								locationsToRemove.add(location);
						}
						
						locationsToAddA.removeAll(locationsToRemove);
					}
					
					
					//second: if we found more than one candidate and at least one of them is overlapping a known signal/splice site, we discard all others
					if(locationsToAddA.size() > 1 && (foundLocationWithKnownSpliceSignal || foundLocationOverlappingKnownJunction)) {
						locationsToRemove.clear();
						for(ReadLocation location : locationsToAddA) {
							if(foundLocationOverlappingKnownJunction) {
								if(!location.overlapsKnownJunction())
									locationsToRemove.add(location);
							}
							else if(foundLocationWithKnownSpliceSignal) {
								if(!location.hasSpliceSignal())
									locationsToRemove.add(location);
							}
						}
						locationsToAddA.removeAll(locationsToRemove);
					}
					
					//third: if there is only one location left we accept it if it overlaps at least with 2bps at each site or has fewer mismatches than the original alignment or overlaps known signal/splice site
					if(locationsToAddA.size() == 1) {
						if(locationsToAddA.get(0).hasSpliceSignal() || locationsToAddA.get(0).overlapsKnownJunction() || locationsToAddA.get(0).getMismatches() < minMismatches)
							locationsToAdd.add(locationsToAddA.get(0));
						else if(locationsToAddA.get(0).getEndA() - locationsToAddA.get(0).getStartA() + 1 > 1 && 
								locationsToAddA.get(0).getEndA() - locationsToAddA.get(0).getStartA() + 1 < (readSequence.length() - 1)) {
							locationsToAdd.add(locationsToAddA.get(0));
						}
					}
					
					//otherwise: the overlap must be larger at each site....
					else {
						for(ReadLocation location : locationsToAddA) {
							if(location.getMismatches() < minMismatches) {
								locationsToAdd.add(location);
							}
							
							else if(location.hasSpliceSignal() || location.overlapsKnownJunction()) {
								if(location.getEndA() - location.getStartA() + 1 > 2 && 
								location.getEndA() - location.getStartA() + 1 < (readSequence.length() - 2)) {
									locationsToAdd.add(location);
								}
							}
							
							else if(location.getEndA() - location.getStartA() + 1 > 4 && 
								location.getEndA() - location.getStartA() + 1 < (readSequence.length() - 4)) {
								locationsToAdd.add(location);
							}
						}
					}
					
					
					prevFullLocationsToAdd.clear();
					prevFullId = currentFullId;
					for(ReadLocation location : locationsToAdd) {
						if(currentAnnotatedSpliceSiteMap == null || !foundLocationOverlappingKnownJunction) {
							if(!this.preferExtensionsWithKnownSpliceSignal || !foundLocationWithKnownSpliceSignal || location.hasSpliceSignal()) {
								//bw.write(String.format("%s\t%s\t%s\t%s\t%s,%s,%s,%s\t%s\t%s\t%s\t%s\t%s\t%s\n", context.getId(),splittedLine[1],chr,strand,location.getStartA(),location.getEndA(),location.getStartB(),location.getEndB(),location.getMismatches(),location.getStrandOfSpliceSignal(),location.overlapsKnownJunction(),splittedLine[11],splittedLine[12],splittedLine[13]));
								sb.setLength(0);
								sb.append(context.getId()).append("\t").append(splittedLine[1]).append("\t").append(chr).append("\t").append(strand).append("\t").append(location.getStartA()).append(",").append(location.getEndA()).append(",").append(location.getStartB()).append(",").append(location.getEndB()).append("\t").append(location.getMismatches()).append("\t").append(location.getStrandOfSpliceSignal()).append("\t").append(location.overlapsKnownJunction()).append("\t").append(splittedLine[11]).append("\t").append(splittedLine[12]).append("\t").append(splittedLine[13]).append("\n");
								bw.write(sb.toString());
								
								ReadLocation l = new ReadLocation(chr,strand,'F',location.getStartA(),location.getEndA(),location.getStartB(),location.getEndB(),location.getMismatches());
								l.setStrandOfSpliceSignal(location.getStrandOfSpliceSignal());
								l.setOverlapsKnownJunction(location.overlapsKnownJunction());
								prevFullLocationsToAdd.add(l);
								
								
							}
						}
						else if(location.overlapsKnownJunction()) {
							sb.setLength(0);
							sb.append(context.getId()).append("\t").append(splittedLine[1]).append("\t").append(chr).append("\t").append(strand).append("\t").append(location.getStartA()).append(",").append(location.getEndA()).append(",").append(location.getStartB()).append(",").append(location.getEndB()).append("\t").append(location.getMismatches()).append("\t").append(location.getStrandOfSpliceSignal()).append("\t").append(location.overlapsKnownJunction()).append("\t").append(splittedLine[11]).append("\t").append(splittedLine[12]).append("\t").append(splittedLine[13]).append("\n");
							bw.write(sb.toString());
							
							ReadLocation l = new ReadLocation(chr,strand,'F',location.getStartA(),location.getEndA(),location.getStartB(),location.getEndB(),location.getMismatches());
							l.setStrandOfSpliceSignal(location.getStrandOfSpliceSignal());
							l.setOverlapsKnownJunction(location.overlapsKnownJunction());
							prevFullLocationsToAdd.add(l);
						}
					}
					
				}
				
				
				//processing split read
				else {
					
					if(prevSplitId != null && currentSplitId.equals(prevSplitId)) {
						for(String line : prevSplitLocationsToAdd) {
							splittedLocation = tabPattern.split(line);
							splittedLocation[1] = splittedLine[1];
							
							bw.write(splittedLocation[0]);
							for(int i = 1; i < splittedLocation.length;i++) {
								bw.write("\t" + splittedLocation[i]);
							}
							bw.write("\t" + splittedLine[11] + "\t" + splittedLine[12] + "\t" + splittedLine[13] + "\n");
						}
						continue;
					}
					
					prevSplitId = currentSplitId;
					prevSplitLocationsToAdd.clear();
					
					//first check if there could be an indel near the splice sites
					//we only do this if the split does not overlap with a known junction or splice signal and has more than 0 mismatches
					if(startB - endA >= 25 && mismatchesOfOriginalAlignment > 0 && strandOfSpliceSignal == '0' && !overlapsKnownJunction) {
						//get the possible deletion sizes
						currentCombinations.clear();
						tmpMinMismatches = Integer.MAX_VALUE;
						shift = 1;
						//check if a deletion occurs upstream of the split
						for(int i = 1; i <= this.maxMissmatches; i++) { 
							tmpStrandOfSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA + i,startB,strand,context.isStrandSpecific());
							if(tmpStrandOfSpliceSignal == '+' || tmpStrandOfSpliceSignal == '-') {
								int newEndA = endA - (2* this.maxMissmatches) - 1;
								if(newEndA < startA) newEndA = startA;
									
								while(newEndA < endA) {
									ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
									coordinates.add(new Pair<Integer,Integer>(startA,newEndA));
									coordinates.add(new Pair<Integer,Integer>(newEndA + i + 1, endA + i));
									coordinates.add(new Pair<Integer,Integer>(startB,endB));
									
									tmpLocation = new ReadLocation(chr,strand,'S',coordinates,tmpMismatches);
									tmpMismatches = this.getMismatchCount(tmpLocation, readSequence,localReference, contextOffset);
									tmpLocation.setMismatches(tmpMismatches);
									tmpLocation.setStrandOfSpliceSignal(tmpStrandOfSpliceSignal);
									tmpLocation.setOverlapsKnownJunction(false);
									if(tmpMismatches < mismatchesOfOriginalAlignment && tmpMismatches < tmpMinMismatches) {
										
										tmpMinMismatches = tmpMismatches;
										currentCombinations.clear();
										currentCombinations.add(tmpLocation);
									}
									newEndA++;
									
								}
							}
						}
						
						
						//check if a deletion occurs downstream of the split
						for(int i = 1; i <= this.maxMissmatches; i++) { 
							tmpStrandOfSpliceSignal = hasSpliceSignal(localReference,contextOffset,donorSite,acceptorSite,endA,startB - i,strand,context.isStrandSpecific());
							if(tmpStrandOfSpliceSignal == '+' || tmpStrandOfSpliceSignal == '-') {
								int newEndB = startB - i;
									
								while(newEndB < startB + (2 * this.maxMissmatches) && newEndB + i < endB) {
									ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
									coordinates.add(new Pair<Integer,Integer>(startA,endA));
									coordinates.add(new Pair<Integer,Integer>(startB - i, newEndB));
									coordinates.add(new Pair<Integer,Integer>(newEndB + i + 1,endB));
									
									tmpLocation = new ReadLocation(chr,strand,'S',coordinates,tmpMismatches);
									tmpMismatches = this.getMismatchCount(tmpLocation, readSequence,localReference, contextOffset);
									tmpLocation.setMismatches(tmpMismatches);
									tmpLocation.setStrandOfSpliceSignal(tmpStrandOfSpliceSignal);
									tmpLocation.setOverlapsKnownJunction(false);
									if(tmpMismatches < mismatchesOfOriginalAlignment && tmpMismatches < tmpMinMismatches) {
										tmpMinMismatches = tmpMismatches;
										currentCombinations.clear();
										currentCombinations.add(tmpLocation);
									}
									
									newEndB++;
									
									
								}
							}
						}
						
						//found just one valid combination. add it to the output. since this combination overlaps a known splice signal and has fewer mismatches than the old combination,
						//we won't check whether the old combination could be extended to a multi split candidate....
						if(currentCombinations.size() == 1) {
							tmpBuilder.setLength(0);
							tmpBuilder.append(context.getId()).append("\t").append(splittedLine[1]).append("\t").append(chr).append("\t").append(strand);
							tmpBuilder.append("\t").append(currentCombinations.get(0).getCoordinates().get(0).getFirst()).append(",").append(currentCombinations.get(0).getCoordinates().get(0).getSecond());
							for(int i = 1; i < currentCombinations.get(0).getCoordinates().size(); i++) {
								tmpBuilder.append(",").append(currentCombinations.get(0).getCoordinates().get(i).getFirst()).append(",").append(currentCombinations.get(0).getCoordinates().get(i).getSecond());
							}
							
							tmpBuilder.append("\t").append(currentCombinations.get(0).getMismatches()).append("\t").append(currentCombinations.get(0).getStrandOfSpliceSignal()).append("\t").append(currentCombinations.get(0).overlapsKnownJunction());
							bw.write(tmpBuilder.toString() + "\t" + splittedLine[11] + "\t" + splittedLine[12] + "\t" + splittedLine[13] + "\n");
							
							prevSplitLocationsToAdd.add(tmpBuilder.toString());
							continue;
						}
					
					}
					
					//no indel found, check if the split can be extended to a multi splti candidate
					if(endA < startB) {
						locationsToAddA = getSplitLocations(startA, endA, readSequence.substring(0,endA-startA + 1),localReference,mismatchesArrayForDownstreamSplits,mismatchesArrayForUpstreamSplits,contextOffset,chr, strand, hashedSplits,insertionKeys,currentAnnotatedSpliceSiteMap,readLocationPool,mismatchArrayPool, context.isStrandSpecific(),true,false);
						
						readLocationPool.clear();
						locationsToAddB = getSplitLocations(startB, endB, readSequence.substring((endA-startA + 1)),localReference,mismatchesArrayForDownstreamSplits,mismatchesArrayForUpstreamSplits,contextOffset,chr, strand, hashedSplits,insertionKeys,currentAnnotatedSpliceSiteMap,readLocationPool,mismatchArrayPool, context.isStrandSpecific(),false,true);
						
						
					}
					
					else {
						if(startB < startA || endB == endA)
							continue;
						
						locationsToAddA = getSplitLocations(startA, startB, readSequence.substring(0,startB-startA + 1),localReference,mismatchesArrayForDownstreamSplits,mismatchesArrayForUpstreamSplits,contextOffset,chr, strand, hashedSplits,insertionKeys,currentAnnotatedSpliceSiteMap,readLocationPool,mismatchArrayPool, context.isStrandSpecific(),true,false);
						
						readLocationPool.clear();
						locationsToAddB = getSplitLocations(endA + 1, endB, readSequence.substring((endA-startA + 1)).substring(endA - startB + 1),localReference,mismatchesArrayForDownstreamSplits,mismatchesArrayForUpstreamSplits,contextOffset,chr, strand, hashedSplits,insertionKeys,currentAnnotatedSpliceSiteMap,readLocationPool, mismatchArrayPool, context.isStrandSpecific(),false,true);
					}
					
					currentCombinations.clear();
					
					for(ReadLocation l1 : locationsToAddA) {
						if(l1.getStartA() < startA + this.maxMissmatches) {
							
							ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
							coordinates.add(new Pair<Integer,Integer>(l1.getStartA(),l1.getEndA()));

							if(endA < startB) {
								coordinates.add(new Pair<Integer,Integer>(l1.getStartB(),l1.getEndB()));
							}
							else {
								if(l1.getEndB() == startB) {
									coordinates.add(new Pair<Integer,Integer>(l1.getStartB(),endA));
								}
								else {
									continue;
								}
							}
							coordinates.add(new Pair<Integer,Integer>(startB,endB));
							
							
							//check if we've created an artificial split
							artificialSplit = false;
							for(int i = 1; i < coordinates.size(); i++) {
								if(coordinates.get(i).getFirst() - 1 == coordinates.get(i-1).getSecond()) {
									artificialSplit = true;
									break;
								}
							}
							if(artificialSplit) {
								continue;
							}
							
							
							tmpLocation = new ReadLocation(chr,strand,'S',coordinates,tmpMismatches);
							
							tmpLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
							tmpLocation.setOverlapsKnownJunction(overlapsKnownJunction);
							
							
							tmpMismatches = this.getMismatchCount(tmpLocation, readSequence,localReference, contextOffset);
							
							tmpLocation.setMismatches(tmpMismatches);
							if(tmpMismatches <= minMismatches + this.maxMissmatchDifference && tmpMismatches <= this.maxMissmatches) {
								if(tmpMismatches < mismatchesOfOriginalAlignment || (l1.getEndA() - l1.getStartA() + 1 > 1)) {
									
									if(tmpMismatches < minMismatches + this.maxMissmatchDifference) {
										currentCombinations.clear();
									}
									
									if(tmpMismatches < minMismatches)
										minMismatches = tmpMismatches;
									
									currentCombinations.add(tmpLocation);
								}
							}
						}
					}
					
				 
					for(ReadLocation l2 : locationsToAddB) {
						if(l2.getStartB() > endB - this.maxMissmatches) {
							ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
							coordinates.add(new Pair<Integer,Integer>(startA,endA));
							if(endA < startB)
								coordinates.add(new Pair<Integer,Integer>(l2.getStartA(),l2.getEndA()));
							else {
								//coordinates.add(new Pair<Integer,Integer>(l2.getStartA() - 1,l2.getEndA()));
								if(l2.getStartA() == endA + 1)
									coordinates.add(new Pair<Integer,Integer>(startB,l2.getEndA()));
								else
									continue;
							}
							
							coordinates.add(new Pair<Integer,Integer>(l2.getStartB(),l2.getEndB()));
							
						
							//check if we've created an artificial split
							artificialSplit = false;
							for(int i = 1; i < coordinates.size(); i++) {
								if(coordinates.get(i).getFirst() - 1 == coordinates.get(i-1).getSecond()) {
									artificialSplit = true;
									break;
								}
							}
							if(artificialSplit) {
								continue;
							}
							
							
							tmpLocation = new ReadLocation(chr,strand,'S',coordinates,tmpMismatches);
							
							tmpLocation.setStrandOfSpliceSignal(strandOfSpliceSignal);
							tmpLocation.setOverlapsKnownJunction(overlapsKnownJunction);
							
							
							tmpMismatches = this.getMismatchCount(tmpLocation, readSequence,localReference, contextOffset);
							tmpLocation.setMismatches(tmpMismatches);
							if(tmpMismatches <= minMismatches + this.maxMissmatchDifference && tmpMismatches <= this.maxMissmatches) {
								if(tmpMismatches < mismatchesOfOriginalAlignment || (l2.getEndB() - l2.getStartB() + 1 > 1)) {
									if(tmpMismatches < minMismatches + this.maxMissmatchDifference) {
										currentCombinations.clear();
									}
									
									if(tmpMismatches < minMismatches)
										minMismatches = tmpMismatches;
									
									currentCombinations.add(tmpLocation);
								}
								
							}
						}
					}
				
					foundLocationWithKnownSpliceSignal = false;
					foundLocationOverlappingKnownJunction = false;
					
					for(ReadLocation location : currentCombinations) {
						if(location.hasSpliceSignal())
							foundLocationWithKnownSpliceSignal = true;
						if(location.overlapsKnownJunction())
							foundLocationOverlappingKnownJunction = true;
					}
					
					for(ReadLocation location : currentCombinations) {
						
						if(location.getMismatches() - minMismatches <= this.maxMissmatchDifference) {
							if(currentAnnotatedSpliceSiteMap == null || !foundLocationOverlappingKnownJunction) {
								if(!this.preferExtensionsWithKnownSpliceSignal || !foundLocationWithKnownSpliceSignal || location.hasSpliceSignal()) {
									tmpBuilder.setLength(0);
									tmpBuilder.append(context.getId()).append("\t").append(splittedLine[1]).append("\t").append(chr).append("\t").append(strand);
									tmpBuilder.append("\t").append(location.getCoordinates().get(0).getFirst()).append(",").append(location.getCoordinates().get(0).getSecond());
									for(int i = 1; i < location.getCoordinates().size(); i++) {
										tmpBuilder.append(",").append(location.getCoordinates().get(i).getFirst()).append(",").append(location.getCoordinates().get(i).getSecond());
									}
									tmpBuilder.append("\t").append(location.getMismatches()).append("\t").append(location.getStrandOfSpliceSignal()).append("\t").append(location.overlapsKnownJunction());
									bw.write(tmpBuilder.toString() + "\t" + splittedLine[11] + "\t" + splittedLine[12] + "\t" + splittedLine[13] + "\n");
									prevSplitLocationsToAdd.add(tmpBuilder.toString());
									
								}
							}
							else if(location.overlapsKnownJunction()) {
								tmpBuilder.setLength(0);
								tmpBuilder.append(context.getId()).append("\t").append(splittedLine[1]).append("\t").append(chr).append("\t").append(strand);
								tmpBuilder.append("\t").append(location.getCoordinates().get(0).getFirst()).append(",").append(location.getCoordinates().get(0).getSecond());
								for(int i = 1; i < location.getCoordinates().size(); i++) {
									tmpBuilder.append(",").append(location.getCoordinates().get(i).getFirst()).append(",").append(location.getCoordinates().get(i).getSecond());
								}
								tmpBuilder.append("\t").append(location.getMismatches()).append("\t").append(location.getStrandOfSpliceSignal()).append("\t").append(location.overlapsKnownJunction());
								bw.write(tmpBuilder.toString() + "\t" + splittedLine[11] + "\t" + splittedLine[12] + "\t" + splittedLine[13] + "\n");
								prevSplitLocationsToAdd.add(tmpBuilder.toString());
							}
						}
					}
				}
			}
			br.close();
			bw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	
	private void getAlreadyCheckedSplits(InitialRead read,HashSet<String> alreadyCheckedSplits,StringBuilder splitKey) {
		alreadyCheckedSplits.clear();
		for(InitialReadLocation readLocation : read.getLocations()) {
			if(readLocation.getMappingType() == 'S') {
				splitKey.setLength(0);
				splitKey.append(readLocation.getStartA());
				splitKey.append('_');
				splitKey.append(readLocation.getEndA());
				alreadyCheckedSplits.add(splitKey.toString());
			}
		}
	}
	
	
	private ArrayList<int[]> getAlreadyCheckedSplits(InitialRead read, ArrayList<int[]> alreadyCheckedSplits) {
		
		alreadyCheckedSplits.clear();
		
		for(InitialReadLocation readLocation : read.getLocations()) {
			if(readLocation.getMappingType() == 'S') {
				int[] tmpSplit = {readLocation.getStartA(),readLocation.getEndA()};
				alreadyCheckedSplits.add(tmpSplit);
			}
		}
		
		return alreadyCheckedSplits;
	}
	
	
	
	private class Split implements Comparable {
		private int splitPositionA;
		private int splitPositionB;
		private double score;
		
		private boolean hasSpliceSignal;
		private boolean overlapsKnownJunction;
		private boolean isIndel;
		
		private HashSet<Integer> startPositions;
		
		
		public Split(int splitPositionA, int splitPositionB) {
			this.splitPositionA = splitPositionA;
			this.splitPositionB = splitPositionB;
			this.startPositions = new HashSet<Integer>();
		}
		
		
		
		public Split(int splitPositionA, int splitPositionB,boolean hasSpliceSignal,boolean overlapsKnownJunction,boolean isIndel) {
			this.splitPositionA = splitPositionA;
			this.splitPositionB = splitPositionB;
			this.startPositions = new HashSet<Integer>();
			
			this.hasSpliceSignal = hasSpliceSignal;
			this.overlapsKnownJunction = overlapsKnownJunction;
			this.isIndel = isIndel;
			
			
		}
		
		
		public boolean hasSpliceSignal() {
			return this.hasSpliceSignal;
		}
		
		public boolean overlapsKnownJunction() {
			return this.overlapsKnownJunction;
		}
		
		public boolean isIndel() {
			return this.isIndel;
		}
		
		public double getScore() {
			return this.score;
		}
		
		public void setScore(double score) {
			this.score = score;
		}
		
		public int getSplitPositionA() {
			return this.splitPositionA;
		}
		
		public int getSplitPositionB() {
			return this.splitPositionB;
		}
		
		public void addStartPosition(int pos) {
			this.startPositions.add(pos);
		}
		
		public int getDifferentStartPositions() {
			return this.startPositions.size();
		}
		
		
		@Override
		public int compareTo(Object splitCompare) {
			Split splitB = (Split)splitCompare;
			return new Integer(this.splitPositionB).compareTo(splitB.getSplitPositionB());
		}
	}
	
	private class SplitScoreComparator implements Comparator<Split> {

		@Override
		public int compare(Split o1, Split o2) {
			// TODO Auto-generated method stub
			return new Double(o2.getScore()).compareTo(o1.getScore());
		}
		
	}
	
	
	private void addAnnotatedSpliceSitesToHash(Context context, TreeMap<Integer,HashSet<Integer>> hashedSplits) {
		//in case we have given an annotation we also hash the predefined junctions here
		HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentMap = null;
		if(!context.isStrandSpecific())
			currentMap = this.chr2annotatedSpliceSites;
		if(context.isStrandSpecific() && context.getStrand().equals("+"))
			currentMap = this.chr2annotatedForwardSpliceSites;
		else if(context.isStrandSpecific() && context.getStrand().equals("-"))
			currentMap = this.chr2annotatedReverseSpliceSites;
		
		if(currentMap != null && currentMap.containsKey(context.getChr())) {
			NavigableMap<Integer,HashSet<Integer>> subMap = currentMap.get(context.getChr()).subMap(context.getStart(), true, context.getEnd(), true);
			for(int i : subMap.keySet()) {
				if(hashedSplits.containsKey(i))
					hashedSplits.get(i).addAll(subMap.get(i));
				else 
					hashedSplits.put(i,subMap.get(i));
			}
		}
	}
	
	
		
	public void extractBestMatchingAlignmentsInLocalResolution(Context context,boolean localContext) {
		try {
			//since we are seeking quite often into the file, we do not need a large buffer here.
			BufferedRandomAccessFile upperReader = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
			BufferedRandomAccessFile lowerReader = new BufferedRandomAccessFile(new File(this.multiMappingFilePath),"r",10240);
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(this.outputPath)),10240);
			long filePointerToFirstOccurence = upperReader.getFilePointer();
			String currentLine = upperReader.getNextLine();
			long prevFilePointer = upperReader.getFilePointer();
			long resetFilePointer = prevFilePointer;
			
			if(currentLine == null) {
				Date date = new Date();
				//System.err.println(String.format("[%s]\tWARNING, found local context without any full reads. context id: %s candidates: %s",date.toLocaleString(), context.getId(),context.getContainedReads()));
				upperReader.close();
				lowerReader.close();
				pw.close();
				return;
			}
			
			StringTokenizer st = new StringTokenizer(currentLine,"\t");
			String currentContextId = st.nextToken();
			String currentReadId = st.nextToken();
			String prevReadId = currentReadId;
			HashSet<String> alreadyProcessedReads = new HashSet<String>();
			//chr strand startA endA startB endB
			st.nextToken();
			st.nextToken();
			st.nextToken();
			st.nextToken();
			int startB = Integer.valueOf(st.nextToken());
			boolean isSplitLocation = false;
			if(startB > 0)
				isSplitLocation = true;
			st.nextToken();
			int missmatches = Integer.valueOf(st.nextToken());
			int minMissmatches = missmatches;
			char strandOfSpliceSignal = st.nextToken().charAt(0);
			boolean overlapsKnownJunction = Boolean.valueOf(st.nextToken());
			boolean processedPartialAndFullReads = false;
			boolean reachedEndOfPartialAndFullBlock = false;
			int processedReads = 0;
			
			ArrayList<MultiMappingLine> currentMappings = new ArrayList<MultiMappingLine>();
			ArrayList<MultiMappingLine> mappingLinePool = new ArrayList<MultiMappingLine>();
			MultiMappingLine tmpLine = new MultiMappingLine(currentLine,missmatches,-1,-1,isSplitLocation,strandOfSpliceSignal,overlapsKnownJunction,false);
			currentMappings.add(tmpLine);
			mappingLinePool.add(tmpLine);
			
			while(true) {
				
				if(!processedPartialAndFullReads || reachedEndOfPartialAndFullBlock)
					currentLine = upperReader.getNextLine();
				else
					currentLine = lowerReader.getNextLine();
				
				//in case we reached the end and we reached it NOT with a seek we break here
				if(currentLine == null && !processedPartialAndFullReads) {
					break;
				}
				
				//in case we are currently seeking for partial reads we go back to the previous position
				else if(currentLine == null) {
					upperReader.seek(resetFilePointer);
					reachedEndOfPartialAndFullBlock = true;
					continue;
				}
				
				
				
				st = new StringTokenizer(currentLine,"\t");
				currentContextId = st.nextToken();
				currentReadId = st.nextToken();
				
				//chr strand startA endA startB endB
				st.nextToken();
				st.nextToken();
				st.nextToken();
				st.nextToken();
				startB = Integer.valueOf(st.nextToken());
				isSplitLocation = false;
				if(startB > 0)
					isSplitLocation = true;
				st.nextToken();
				missmatches = Integer.valueOf(st.nextToken());
				strandOfSpliceSignal = st.nextToken().charAt(0);
				overlapsKnownJunction = Boolean.valueOf(st.nextToken());
				
				if(prevReadId.equals(currentReadId)) {
					
					if(currentMappings.size() < mappingLinePool.size()) {
						tmpLine = mappingLinePool.get(currentMappings.size());
						tmpLine.update(currentLine,missmatches,-1,-1,isSplitLocation,strandOfSpliceSignal,overlapsKnownJunction,false);
					}
					else {
						tmpLine = new MultiMappingLine(currentLine,missmatches,-1,-1,isSplitLocation,strandOfSpliceSignal,overlapsKnownJunction,false);
						mappingLinePool.add(tmpLine);
					}
					
					currentMappings.add(tmpLine);
					if(missmatches < minMissmatches)
						minMissmatches = missmatches;
				}
				
				else if(!processedPartialAndFullReads && context.getPartialAndFullReads2filePointer().containsKey(prevReadId) && context.getPartialAndFullReads2filePointer().get(prevReadId) != filePointerToFirstOccurence) {
					processedPartialAndFullReads = true;
					lowerReader.seek(context.getPartialAndFullReads2filePointer().get(prevReadId));
					resetFilePointer = prevFilePointer;
				}
				
				else if(processedPartialAndFullReads && !reachedEndOfPartialAndFullBlock && context.getPartialAndFullReads2filePointer().containsKey(prevReadId)) {
					reachedEndOfPartialAndFullBlock = true;
					//just go back one line with the upper reader (should be buffered)
					upperReader.seek(resetFilePointer);
					
				}
				
				else {
					if(!alreadyProcessedReads.contains(prevReadId)) {
						Collections.sort(currentMappings);
						printRead(currentMappings,minMissmatches,this.maxMissmatchDifference,pw);
						if(context.getPartialAndFullReads2filePointer().containsKey(prevReadId)) {
							context.getPartialAndFullReads2filePointer().remove(prevReadId);
							alreadyProcessedReads.add(prevReadId);
						}
					}
					
					minMissmatches = missmatches;
					prevReadId = currentReadId;
					currentMappings.clear();
					
					if(currentMappings.size() < mappingLinePool.size()) {
						tmpLine = mappingLinePool.get(currentMappings.size());
						tmpLine.update(currentLine,missmatches,-1,-1,isSplitLocation,strandOfSpliceSignal,overlapsKnownJunction,false);
					}
					else {
						tmpLine = new MultiMappingLine(currentLine,missmatches,-1,-1,isSplitLocation,strandOfSpliceSignal,overlapsKnownJunction,false);
						mappingLinePool.add(tmpLine);
					}
					currentMappings.add(tmpLine);
					
					processedPartialAndFullReads = false;
					reachedEndOfPartialAndFullBlock = false;
					filePointerToFirstOccurence = prevFilePointer;
					processedReads++;
					
				}
				prevFilePointer = upperReader.getFilePointer();
			}
			//print the last read
			Collections.sort(currentMappings);
			printRead(currentMappings,minMissmatches,this.maxMissmatchDifference,pw);
			
			upperReader.close();
			lowerReader.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	public void extractBestMatchingPairedEndAlignments(int maxContextSize, boolean predictPolyA) {
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(this.multiMappingFilePath)));
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(this.outputPath)));
			String currentLine = br.readLine();
			
			if(currentLine == null) {
				Date date = new Date();
				//System.err.println(String.format("[%s]\tWARNING, found local context without any full reads.",date.toLocaleString()));
				br.close();
				pw.close();
				return;
			}
			
			StringBuffer sb = new StringBuffer();
			StringTokenizer st = new StringTokenizer(currentLine,"\t");
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedCoordinates;
			String contextId = st.nextToken();
			String readId = st.nextToken();
			String readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
			String prevReadId = readId;
			String prevReadIdPrefix = readIdPrefix;
			String chr = st.nextToken();
			char strand = st.nextToken().toCharArray()[0];
			splittedCoordinates = commaPattern.split(st.nextToken());
			ArrayList<Pair<Integer,Integer>> coordinates = new ArrayList<Pair<Integer,Integer>>();
			for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
				coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
			}
			char mappingType;
			if(coordinates.size() == 1 || coordinates.get(1).getFirst() <= 0)
				mappingType = 'F';
			else
				mappingType = 'S';
			
			int mismatches = Integer.valueOf(st.nextToken());
			int minMismatches;
			char hasSpliceSignal = st.nextToken().charAt(0);
			boolean overlapsKnownJunction = Boolean.valueOf(st.nextToken());
			
			double readScore = Double.valueOf(st.nextToken());
			double locationScore = Double.valueOf(st.nextToken());
			
			int overallMappingCount = Integer.valueOf(st.nextToken());
			int overallValidPairCount = Integer.valueOf(st.nextToken());
			char overlapsPolyAtail = st.nextToken().charAt(0);
			
			boolean foundReadWithSpliceSignal;
			boolean foundReadOverlappingKnownJunction;
			boolean removeClippedAlignments;
			
			ArrayList<ReadLocation> currentReadLocations = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> locationsToRemove = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> tmpReadLocations;
			ArrayList<ReadLocation> tmpReadLocationsToPrintA = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> tmpReadLocationsToPrintB = new ArrayList<ReadLocation>();
			ReadLocation tmpReadLocationA;
			ReadLocation tmpReadLocationB;
			ReadLocation currentReadLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
			currentReadLocation.setReadId(readId);
			currentReadLocation.setStrandOfSpliceSignal(hasSpliceSignal);
			currentReadLocation.setOverlapsKnownJunction(overlapsKnownJunction);
			currentReadLocation.setScore(locationScore);
			currentReadLocation.setOverlapsPolyAtail(overlapsPolyAtail);
			currentReadLocations.add(currentReadLocation);
			
			
			Read read = new Read(readId);
			read.setScore(readScore);
			read.setOverallMappingCount(overallMappingCount);
			read.setOverallValidPairCount(overallValidPairCount);
			ReadPair<Read,Read> currentPair = new ReadPair<Read,Read>(readIdPrefix);
			ReadPair<Read,Read> tmpPair = null;
			currentPair.setFirst(read);
			
			while((currentLine = br.readLine()) != null) {
				st = new StringTokenizer(currentLine,"\t");
				contextId = st.nextToken();
				readId = st.nextToken();
				readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
				chr = st.nextToken();
				strand = st.nextToken().toCharArray()[0];
				splittedCoordinates = commaPattern.split(st.nextToken());
				coordinates = new ArrayList<Pair<Integer,Integer>>();
				for(int i = 0; i < splittedCoordinates.length - 1; i+=2) {
					coordinates.add(new Pair<Integer,Integer>(Integer.valueOf(splittedCoordinates[i]),Integer.valueOf(splittedCoordinates[i+1])));
				}
				if(coordinates.size() == 1 || coordinates.get(1).getFirst() <= 0)
					mappingType = 'F';
				else
					mappingType = 'S';

				mismatches = Integer.valueOf(st.nextToken());
				hasSpliceSignal = st.nextToken().charAt(0);
				overlapsKnownJunction = Boolean.valueOf(st.nextToken());
				readScore = Double.valueOf(st.nextToken());
				locationScore = Double.valueOf(st.nextToken());
				overallMappingCount = Integer.valueOf(st.nextToken());
				overallValidPairCount = Integer.valueOf(st.nextToken());
				overlapsPolyAtail = st.nextToken().charAt(0);
				
				if(readId.equals(prevReadId)) {
					currentReadLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
					currentReadLocation.setReadId(readId);
					currentReadLocation.setStrandOfSpliceSignal(hasSpliceSignal);
					currentReadLocation.setOverlapsKnownJunction(overlapsKnownJunction);
					currentReadLocation.setScore(locationScore);
					currentReadLocation.setOverlapsPolyAtail(overlapsPolyAtail);
					currentReadLocations.add(currentReadLocation);
				}
				
				//second read from the actual fragment
				else if(readIdPrefix.equals(prevReadIdPrefix)) {
					//the actual sparse read is readA of the current pair
					for(int i = 0; i < currentReadLocations.size(); i++) {
						read.addLocation(currentReadLocations.get(i));
					}
					
					read = new Read(readId);
					read.setScore(readScore);
					read.setOverallMappingCount(overallMappingCount);
					read.setOverallValidPairCount(overallValidPairCount);
					currentPair.setSecond(read);
					currentReadLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
					currentReadLocation.setReadId(readId);
					currentReadLocation.setStrandOfSpliceSignal(hasSpliceSignal);
					currentReadLocation.setOverlapsKnownJunction(overlapsKnownJunction);
					currentReadLocation.setScore(locationScore);
					currentReadLocation.setOverlapsPolyAtail(overlapsPolyAtail);
					currentReadLocations.clear();
					currentReadLocations.add(currentReadLocation);
					
					prevReadId = readId;
					
				}
				
				//here we have a new fragment
				else {
					
					//first process the old fragment/read
					for(int i = 0; i < currentReadLocations.size(); i++) {
						read.addLocation(currentReadLocations.get(i));
					}
					
					if(currentPair.getSecond() != null) {
						currentPair.getSecond().setLocations(removeArtificialSplitLocations(currentPair.getSecond().getLocations()));
						if(currentPair.getSecond().getLocations().isEmpty())
							currentPair.setSecond(null);
					}
					
					currentPair.getFirst().setLocations(removeArtificialSplitLocations(currentPair.getFirst().getLocations()));
					if(currentPair.getFirst().getLocations().isEmpty()) {
						if(currentPair.getSecond() != null) {
							currentPair.setFirst(currentPair.getSecond());
							currentPair.setSecond(null);
						}
					}
					
					//check for each mate: if we find a full/split alignment we will remove all clipped alignments
					//first mate
					removeClippedAlignments(currentPair.getFirst(),locationsToRemove);
					//second mate
					if(currentPair.getSecond() != null)
						removeClippedAlignments(currentPair.getSecond(),locationsToRemove);
					
					
					//both mates found. check for valid pairs
					if(currentPair.getSecond() != null && !currentPair.getFirst().getLocations().isEmpty()) {
						getValidReadPairs(currentPair,maxContextSize);
											
						if(!currentPair.getValidPairs().isEmpty()) {
							//get best matching mappings from valid pairs
							foundReadWithSpliceSignal = false;
							foundReadOverlappingKnownJunction = false;
							minMismatches = Integer.MAX_VALUE;
							ArrayList<Pair<Integer,Integer>> filteredPairs = new ArrayList<Pair<Integer,Integer>>();
							for(Pair<Integer,Integer> indexPair : currentPair.getValidPairs()) {
								tmpReadLocationA = currentPair.getFirst().getLocations().get(indexPair.getFirst());
								tmpReadLocationB = currentPair.getSecond().getLocations().get(indexPair.getSecond());
								if(tmpReadLocationA.getMismatches() + tmpReadLocationB.getMismatches() < minMismatches) {
									minMismatches = tmpReadLocationA.getMismatches() + tmpReadLocationB.getMismatches();
									filteredPairs.clear();
									filteredPairs.add(indexPair);
									foundReadWithSpliceSignal = false;
									foundReadOverlappingKnownJunction = false;
									
									if(tmpReadLocationA.getMappingType() == 'S') {
										if(tmpReadLocationA.hasSpliceSignal()) foundReadWithSpliceSignal = true;
										if(tmpReadLocationA.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
									}
									
									if(tmpReadLocationB.getMappingType() == 'S') {
										if(tmpReadLocationB.hasSpliceSignal()) foundReadWithSpliceSignal = true;
										if(tmpReadLocationB.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
									}
									
									
									
								}
								
								else if(tmpReadLocationA.getMismatches() + tmpReadLocationB.getMismatches() <= minMismatches + this.maxMissmatchDifference) {
									filteredPairs.add(indexPair);
									if(tmpReadLocationA.getMappingType() == 'S') {
										if(tmpReadLocationA.hasSpliceSignal()) foundReadWithSpliceSignal = true;
										if(tmpReadLocationA.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
									}
									
									if(tmpReadLocationB.getMappingType() == 'S') {
										if(tmpReadLocationB.hasSpliceSignal()) foundReadWithSpliceSignal = true;
										if(tmpReadLocationB.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
									}
								}
							}
							
							
							//print out the valid pairs
							tmpReadLocationsToPrintA.clear();
							tmpReadLocationsToPrintB.clear();
							for(Pair<Integer,Integer> indexPair : filteredPairs) {
								tmpReadLocationA = currentPair.getFirst().getLocations().get(indexPair.getFirst());
								tmpReadLocationB = currentPair.getSecond().getLocations().get(indexPair.getSecond());
								
								if(tmpReadLocationA.getMappingType() == 'F' && tmpReadLocationB.getMappingType() == 'F') {
									tmpReadLocationsToPrintA.add(tmpReadLocationA);
									tmpReadLocationsToPrintB.add(tmpReadLocationB);
									continue;
								}
								
								if(!foundReadWithSpliceSignal && !foundReadOverlappingKnownJunction) {
									tmpReadLocationsToPrintA.add(tmpReadLocationA);
									tmpReadLocationsToPrintB.add(tmpReadLocationB);
									continue;
								}
								
								if(foundReadOverlappingKnownJunction) {
									if(tmpReadLocationA.overlapsKnownJunction() || tmpReadLocationB.overlapsKnownJunction()) {
										tmpReadLocationsToPrintA.add(tmpReadLocationA);
										tmpReadLocationsToPrintB.add(tmpReadLocationB);
									}
									continue;
								}
								
								if(foundReadWithSpliceSignal) {
									if(tmpReadLocationA.hasSpliceSignal() || tmpReadLocationB.hasSpliceSignal()) {
										tmpReadLocationsToPrintA.add(tmpReadLocationA);
										tmpReadLocationsToPrintB.add(tmpReadLocationB);
									}
									continue;
								}
							}
							
							for(ReadLocation l : tmpReadLocationsToPrintA) {
								printPairedEndAlignment(pw,sb,l,currentPair.getFirst().getScore(),currentPair.getFirst().getOverallMappingCount(),currentPair.getFirst().getOverallValidPairCount(),contextId);
							}
							
							for(ReadLocation l : tmpReadLocationsToPrintB) {
								printPairedEndAlignment(pw,sb,l,currentPair.getSecond().getScore(),currentPair.getSecond().getOverallMappingCount(),currentPair.getSecond().getOverallValidPairCount(),contextId);
							}
							
						}
						
						//if both mates mapped, but no valid pair found, we check if there are alignments overlapping a cleavage site
						else if(predictPolyA) {
							for(ReadLocation l : currentPair.getFirst().getLocations()) {
								if(l.getOverlapsPolyAtail() == '1') 
									printPairedEndAlignment(pw,sb,l,currentPair.getFirst().getScore(),currentPair.getFirst().getOverallMappingCount(),currentPair.getFirst().getOverallValidPairCount(),contextId);
							}
							
							for(ReadLocation l : currentPair.getSecond().getLocations()) {
								if(l.getOverlapsPolyAtail() == '1') 
									printPairedEndAlignment(pw,sb,l,currentPair.getSecond().getScore(),currentPair.getSecond().getOverallMappingCount(),currentPair.getSecond().getOverallValidPairCount(),contextId);
							}
						}
						
					}

					
					//found only one mate, but with several mapping positions.
					else {
						minMismatches = Integer.MAX_VALUE;
						for(ReadLocation l : currentPair.getFirst().getLocations()) {
							if(l.getMismatches() < minMismatches)
								minMismatches = l.getMismatches();
						}
						tmpReadLocations = filterForMinMismatchesAndKnownSignals(currentPair.getFirst().getLocations(), minMismatches);
						
						for(ReadLocation l : tmpReadLocations) {
							printPairedEndAlignment(pw,sb,l,currentPair.getFirst().getScore(),currentPair.getFirst().getOverallMappingCount(),currentPair.getFirst().getOverallValidPairCount(),contextId);
						}
					}
					
					
					prevReadId = readId;
					prevReadIdPrefix = readIdPrefix;
					
					read = new Read(readId);
					read.setScore(readScore);
					read.setOverallMappingCount(overallMappingCount);
					read.setOverallValidPairCount(overallValidPairCount);
					currentPair = new ReadPair<Read,Read>(readIdPrefix);
					currentPair.setFirst(read);
					
					
					currentReadLocation = new ReadLocation(chr, strand, mappingType, coordinates, mismatches);
					currentReadLocation.setReadId(readId);
					currentReadLocation.setStrandOfSpliceSignal(hasSpliceSignal);
					currentReadLocation.setOverlapsKnownJunction(overlapsKnownJunction);
					currentReadLocation.setScore(locationScore);
					currentReadLocation.setOverlapsPolyAtail(overlapsPolyAtail);
					currentReadLocations.clear();
					currentReadLocations.add(currentReadLocation);
					
				}
			}
			
			
			//proces last bunch here
			for(int i = 0; i < currentReadLocations.size(); i++) {
				read.addLocation(currentReadLocations.get(i));
			}
			
			if(currentPair.getSecond() != null) {
				currentPair.getSecond().setLocations(removeArtificialSplitLocations(currentPair.getSecond().getLocations()));
				if(currentPair.getSecond().getLocations().isEmpty())
					currentPair.setSecond(null);
			}
			
			currentPair.getFirst().setLocations(removeArtificialSplitLocations(currentPair.getFirst().getLocations()));
			if(currentPair.getFirst().getLocations().isEmpty()) {
				if(currentPair.getSecond() != null) {
					currentPair.setFirst(currentPair.getSecond());
					currentPair.setSecond(null);
				}
			}
			
			//check for each mate: if we find a full/split alignment we will remove all clipped alignments
			//first mate
			removeClippedAlignments(currentPair.getFirst(),locationsToRemove);
			//second mate
			if(currentPair.getSecond() != null)
				removeClippedAlignments(currentPair.getSecond(),locationsToRemove);
			
			
			//both mates found. check for valid pairs
			if(currentPair.getSecond() != null && !currentPair.getFirst().getLocations().isEmpty()) {
				
				//if(currentPair.getFirst().getLocations().size() > 1 || currentPair.getSecond().getLocations().size() > 1)
					getValidReadPairs(currentPair,maxContextSize);
				
					if(!currentPair.getValidPairs().isEmpty()) {
					//get best matching mappings from valid pairs
					foundReadWithSpliceSignal = false;
					foundReadOverlappingKnownJunction = false;
					minMismatches = Integer.MAX_VALUE;
					ArrayList<Pair<Integer,Integer>> filteredPairs = new ArrayList<Pair<Integer,Integer>>();
					for(Pair<Integer,Integer> indexPair : currentPair.getValidPairs()) {
						tmpReadLocationA = currentPair.getFirst().getLocations().get(indexPair.getFirst());
						tmpReadLocationB = currentPair.getSecond().getLocations().get(indexPair.getSecond());
						if(tmpReadLocationA.getMismatches() + tmpReadLocationB.getMismatches() < minMismatches) {
							minMismatches = tmpReadLocationA.getMismatches() + tmpReadLocationB.getMismatches();
							filteredPairs.clear();
							filteredPairs.add(indexPair);
							foundReadWithSpliceSignal = false;
							foundReadOverlappingKnownJunction = false;
							
							if(tmpReadLocationA.getMappingType() == 'S') {
								if(tmpReadLocationA.hasSpliceSignal()) foundReadWithSpliceSignal = true;
								if(tmpReadLocationA.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
							}
							
							if(tmpReadLocationB.getMappingType() == 'S') {
								if(tmpReadLocationB.hasSpliceSignal()) foundReadWithSpliceSignal = true;
								if(tmpReadLocationB.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
							}
							
							
							
						}
						
						else if(tmpReadLocationA.getMismatches() + tmpReadLocationB.getMismatches() <= minMismatches + this.maxMissmatchDifference) {
							filteredPairs.add(indexPair);
							if(tmpReadLocationA.getMappingType() == 'S') {
								if(tmpReadLocationA.hasSpliceSignal()) foundReadWithSpliceSignal = true;
								if(tmpReadLocationA.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
							}
							
							if(tmpReadLocationB.getMappingType() == 'S') {
								if(tmpReadLocationB.hasSpliceSignal()) foundReadWithSpliceSignal = true;
								if(tmpReadLocationB.overlapsKnownJunction()) foundReadOverlappingKnownJunction = true;
							}
						}
					}
					
					
					//print out the valid pairs
					tmpReadLocationsToPrintA.clear();
					tmpReadLocationsToPrintB.clear();
					for(Pair<Integer,Integer> indexPair : filteredPairs) {
						tmpReadLocationA = currentPair.getFirst().getLocations().get(indexPair.getFirst());
						tmpReadLocationB = currentPair.getSecond().getLocations().get(indexPair.getSecond());
						
						if(tmpReadLocationA.getMappingType() == 'F' && tmpReadLocationB.getMappingType() == 'F') {
							tmpReadLocationsToPrintA.add(tmpReadLocationA);
							tmpReadLocationsToPrintB.add(tmpReadLocationB);
							continue;
						}
						
						if(!foundReadWithSpliceSignal && !foundReadOverlappingKnownJunction) {
							tmpReadLocationsToPrintA.add(tmpReadLocationA);
							tmpReadLocationsToPrintB.add(tmpReadLocationB);
							continue;
						}
						
						if(foundReadOverlappingKnownJunction) {
							if(tmpReadLocationA.overlapsKnownJunction() || tmpReadLocationB.overlapsKnownJunction()) {
								tmpReadLocationsToPrintA.add(tmpReadLocationA);
								tmpReadLocationsToPrintB.add(tmpReadLocationB);
							}
							continue;
						}
						
						if(foundReadWithSpliceSignal) {
							if(tmpReadLocationA.hasSpliceSignal() || tmpReadLocationB.hasSpliceSignal()) {
								tmpReadLocationsToPrintA.add(tmpReadLocationA);
								tmpReadLocationsToPrintB.add(tmpReadLocationB);
							}
							continue;
						}
					}
					
					for(ReadLocation l : tmpReadLocationsToPrintA) {
						printPairedEndAlignment(pw,sb,l,currentPair.getFirst().getScore(),currentPair.getFirst().getOverallMappingCount(),currentPair.getFirst().getOverallValidPairCount(),contextId);
					}
					
					for(ReadLocation l : tmpReadLocationsToPrintB) {
						printPairedEndAlignment(pw,sb,l,currentPair.getSecond().getScore(),currentPair.getSecond().getOverallMappingCount(),currentPair.getSecond().getOverallValidPairCount(),contextId);
					}
				}
					
				//if both mates mapped, but no valid pair found, we check if one mate overlaps a cleavage site
				else if(predictPolyA) {
					for(ReadLocation l : currentPair.getFirst().getLocations()) {
						if(l.getOverlapsPolyAtail() == '1') 
							printPairedEndAlignment(pw,sb,l,currentPair.getFirst().getScore(),currentPair.getFirst().getOverallMappingCount(),currentPair.getFirst().getOverallValidPairCount(),contextId);
					}
					
					for(ReadLocation l : currentPair.getSecond().getLocations()) {
						if(l.getOverlapsPolyAtail() == '1') 
							printPairedEndAlignment(pw,sb,l,currentPair.getSecond().getScore(),currentPair.getSecond().getOverallMappingCount(),currentPair.getSecond().getOverallValidPairCount(),contextId);
					}
				}
				
			}

			
			//found only one mate, but with several mapping positions.
			else {
				minMismatches = Integer.MAX_VALUE;
				for(ReadLocation l : currentPair.getFirst().getLocations()) {
					if(l.getMismatches() < minMismatches)
						minMismatches = l.getMismatches();
				}
				tmpReadLocations = filterForMinMismatchesAndKnownSignals(currentPair.getFirst().getLocations(), minMismatches);
				
				for(ReadLocation l : tmpReadLocations) {
					printPairedEndAlignment(pw,sb,l,currentPair.getFirst().getScore(),currentPair.getFirst().getOverallMappingCount(),currentPair.getFirst().getOverallValidPairCount(),contextId);
				}
			}
			
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	//removes clipped alignments if we find at least one full or split alignment without clipping
	private void removeClippedAlignments(Read read, ArrayList<ReadLocation> locationsToRemove) {
		locationsToRemove.clear();
		boolean removeClippedAlignments = false;
		for(ReadLocation location : read.getLocations()) {
			if(location.getMappingType() == 'F' && location.getCoordinates().size() > 1 && (location.getCoordinates().get(1).getFirst() < 0 ||location.getCoordinates().get(1).getSecond() < 0))
				locationsToRemove.add(location);
			else
				removeClippedAlignments = true;
		}
		if(removeClippedAlignments) {
			read.getLocations().removeAll(locationsToRemove);
		}
	}
	
	private ArrayList<ReadLocation> removeArtificialSplitLocations(ArrayList<ReadLocation> locations) {
		ArrayList<ReadLocation> filteredLocations = new ArrayList<ReadLocation>();
		ArrayList<ReadLocation> splitLocations = new ArrayList<ReadLocation>();
		boolean hasFullAlignment = false;
		int fullAlignmentMinMismatches = Integer.MAX_VALUE;
		int splitAlignmentMinMismatches = Integer.MAX_VALUE;
		for(ReadLocation l : locations) {
			if(l.getMappingType() == 'F') {
				hasFullAlignment = true;
				if(l.getMismatches() < fullAlignmentMinMismatches)
					fullAlignmentMinMismatches = l.getMismatches();
			}
			
			//only consider single split alignments
			if(l.getMappingType() == 'F' || (l.getCoordinates().size() > 2 && l.getCoordinates().get(l.getCoordinates().size() - 1).getFirst() > 0))
				filteredLocations.add(l);
			else {
				splitLocations.add(l);
				if(l.getMismatches() < splitAlignmentMinMismatches) {
					splitAlignmentMinMismatches = l.getMismatches();
				}
			}
		}
		
		
		if((hasFullAlignment && splitAlignmentMinMismatches >= fullAlignmentMinMismatches && splitLocations.size() > 1) ||
		   (!hasFullAlignment && splitLocations.size() > 3)) {
			HashMap<Integer,MutableInt> splitSize2frequency = new HashMap<Integer,MutableInt>();
			int splitSize;
			for(ReadLocation l : splitLocations) {
				splitSize = l.getCoordinates().get(0).getSecond() - l.getCoordinates().get(0).getFirst() + 1;
				if(splitSize2frequency.containsKey(splitSize))
					splitSize2frequency.get(splitSize).increment();
				else
					splitSize2frequency.put(splitSize, new MutableInt(1));
				
			}
			
			
			for(MutableInt freq : splitSize2frequency.values()) {
				if(freq.intValue() > splitLocations.size()/2) {
					splitLocations.clear();
					break;
				}
			}
		}
		
		filteredLocations.addAll(splitLocations);
		return filteredLocations;
	}
	
	
	
	private void printPairedEndAlignment(UnsynchronizedBufferedWriter pw, StringBuffer sb, ReadLocation l, double readScore, int overallMappingCount, int overallValidPairCount, String contextId) throws Exception {
		sb.setLength(0);
		sb.append(contextId).append("\t").append(l.getReadId()).append("\t").append(l.getChr()).append("\t").append(l.getStrand()).append("\t").append(l.getCoordinates().get(0).getFirst()).append(",").append(l.getCoordinates().get(0).getSecond());
		for(int i = 1; i < l.getCoordinates().size();i++) {
			sb.append(",").append(l.getCoordinates().get(i).getFirst()).append(",").append(l.getCoordinates().get(i).getSecond());
		}
		sb.append("\t").append(l.getMismatches()).append("\t").append(l.getStrandOfSpliceSignal()).append("\t").append(l.overlapsKnownJunction()).append("\t").append(readScore).append("\t").append(l.getScore()).append("\t").append(overallMappingCount).append("\t").append(overallValidPairCount).append("\t").append(l.getOverlapsPolyAtail());
		
		pw.write(sb.toString() + "\n");
	}
	
	
	
	
	private ArrayList<ReadLocation> filterForMinMismatchesAndKnownSignals(ArrayList<ReadLocation> readLocations, int minMismatches) {
		ArrayList<ReadLocation> filteredMappings = new ArrayList<ReadLocation>();
		boolean foundMappingWithKnownSpliceSignal = false;
		boolean foundMappingOverlappingKnownJunction = false;
		
		
		for(int i  = 0;i< readLocations.size(); i++) {
			if(readLocations.get(i).getMismatches() - minMismatches <= this.maxMissmatchDifference) {
				if(readLocations.get(i).hasSpliceSignal()) foundMappingWithKnownSpliceSignal = true;
				if(readLocations.get(i).overlapsKnownJunction()) foundMappingOverlappingKnownJunction = true;
			}
		}
		
		for(int i = 0;i < readLocations.size(); i++) {
			if(readLocations.get(i).getMismatches() - minMismatches <= this.maxMissmatchDifference) {
				if(readLocations.get(i).getMappingType() == 'F') {
					filteredMappings.add(readLocations.get(i));
					continue;
				}
				
				if(!foundMappingOverlappingKnownJunction && (!foundMappingWithKnownSpliceSignal || !this.preferExtensionsWithKnownSpliceSignal)) {
					filteredMappings.add(readLocations.get(i));
					continue;
				}
				
				if(foundMappingOverlappingKnownJunction && readLocations.get(i).overlapsKnownJunction()) {
					filteredMappings.add(readLocations.get(i));
					continue;
				}
				
				if(!foundMappingOverlappingKnownJunction && foundMappingWithKnownSpliceSignal && readLocations.get(i).hasSpliceSignal()) {
					filteredMappings.add(readLocations.get(i));
					continue;
				}
					
			}
		}
		
		return filteredMappings;
	}
	
	
	
	public void extractBestMatchingAlignments(boolean isGlobalResolution) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(this.multiMappingFilePath)));
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(this.outputPath)));
			String currentLine = br.readLine();
			
			if(currentLine == null) {
				Date date = new Date();
				//System.err.println(String.format("[%s]\tWARNING, found local context without any full reads.",date.toLocaleString()));
				br.close();
				pw.close();
				return;
			}
			
			StringTokenizer st = new StringTokenizer(currentLine,"\t");
			String currentContextId = st.nextToken();
			String currentReadId = st.nextToken();
			String prevReadId = currentReadId;
			//chr strand startA endA startB endB
			st.nextToken();
			st.nextToken();
			
			ArrayList<Pair<Integer,Integer>> coordinates;
			Pattern commaPattern = Pattern.compile(",");
			String[] splittedCoordinates = commaPattern.split(st.nextToken());
			boolean isClippedLocation = false;
			boolean isSplitLocation = false;
			if(splittedCoordinates.length > 2) {
				if(Integer.valueOf(splittedCoordinates[splittedCoordinates.length - 2]) > 0)
					isSplitLocation = true;
				
				else if(Integer.valueOf(splittedCoordinates[splittedCoordinates.length - 2]) < 0 || Integer.valueOf(splittedCoordinates[splittedCoordinates.length - 1]) < 0)
					isClippedLocation = true;
				
			}


			int missmatches = Integer.valueOf(st.nextToken());
			char hasKnownSpliceSignal = st.nextToken().charAt(0);
			boolean overlapsKnownJunction = Boolean.valueOf(st.nextToken());
			st.nextToken();
			st.nextToken();
			boolean overlapsCleavageSite = st.nextToken().charAt(0) == '1';
			boolean foundCandidateSplit = isSplitLocation;
			
			
			
			int minMissmatches = Integer.MAX_VALUE;
			if(missmatches < minMissmatches && !isClippedLocation)
				minMissmatches = missmatches;
			
			
			
			int processedReads = 0;
			
			ArrayList<MultiMappingLine> currentMappings = new ArrayList<MultiMappingLine>();
			ArrayList<MultiMappingLine> currentClippedMappings = new ArrayList<MultiMappingLine>();
			ArrayList<MultiMappingLine> mappingLinePool = new ArrayList<MultiMappingLine>();
			MultiMappingLine tmpLine = new MultiMappingLine(currentLine,missmatches,-1,-1,isSplitLocation,hasKnownSpliceSignal,overlapsKnownJunction,overlapsCleavageSite);
			
			if(!isClippedLocation) {
				currentMappings.add(tmpLine);
			}
			else {
				currentClippedMappings.add(tmpLine);
			}
			
			mappingLinePool.add(tmpLine);
			
			while((currentLine = br.readLine()) != null) {
				st = new StringTokenizer(currentLine,"\t");
				currentContextId = st.nextToken();
				currentReadId = st.nextToken();
				//chr strand startA endA startB endB
				st.nextToken();
				st.nextToken();
				splittedCoordinates = commaPattern.split(st.nextToken());
				
				isSplitLocation = false;
				isClippedLocation = false;
				if(splittedCoordinates.length > 2) {
					if(Integer.valueOf(splittedCoordinates[splittedCoordinates.length - 2]) > 0)
						isSplitLocation = true;
					
					else if(Integer.valueOf(splittedCoordinates[splittedCoordinates.length - 2]) < 0 || Integer.valueOf(splittedCoordinates[splittedCoordinates.length - 1]) < 0)
						isClippedLocation = true;
				}
				
				missmatches = Integer.valueOf(st.nextToken());
				hasKnownSpliceSignal = st.nextToken().charAt(0);
				overlapsKnownJunction = Boolean.valueOf(st.nextToken());
				st.nextToken();
				st.nextToken();
				overlapsCleavageSite = st.nextToken().charAt(0) == '1';
				
				if(prevReadId.equals(currentReadId)) {
					
					if(currentMappings.size() + currentClippedMappings.size() < mappingLinePool.size()) {
						tmpLine = mappingLinePool.get(currentMappings.size() + currentClippedMappings.size());
						tmpLine.update(currentLine,missmatches,-1,-1,isSplitLocation,hasKnownSpliceSignal,overlapsKnownJunction,overlapsCleavageSite);
					}
					else {
						tmpLine = new MultiMappingLine(currentLine,missmatches,-1,-1,isSplitLocation,hasKnownSpliceSignal,overlapsKnownJunction,overlapsCleavageSite);
						mappingLinePool.add(tmpLine);
					}
					
					if(!isClippedLocation) {
						currentMappings.add(tmpLine);
					}
					else {
						currentClippedMappings.add(tmpLine);
					}
					
					
					if(missmatches < minMissmatches && !isClippedLocation)
						minMissmatches = missmatches;
					
					if(isSplitLocation)
						foundCandidateSplit = true;
					
				}
			
				else {
					
					if(isGlobalResolution && !currentClippedMappings.isEmpty()) {
						if(currentMappings.isEmpty()) {
							currentMappings.addAll(currentClippedMappings);
							minMissmatches = Integer.MAX_VALUE;
							for(MultiMappingLine ml : currentMappings) {
								if(ml.getMissmatches() < minMissmatches)
									minMissmatches = ml.getMissmatches();
							}
						}
					}
					
					
					else if(!isGlobalResolution && !currentClippedMappings.isEmpty()) {
						//the situation in a default ContextMap run
						int prevSize = currentMappings.size();
						if(currentMappings.isEmpty()) {
							currentMappings.addAll(currentClippedMappings);
						}
						//clipped and regular alignments only appear in polyA tail recognition mode
						//if we found a candidate split alignment we will discard all clipped alignments.
						//otherwise we will add only clipped alignment overlapping with a cleavage site
						else if(!foundCandidateSplit) {
							for(MultiMappingLine ml : currentClippedMappings) {
								if(tmpLine.overlapsCleavageSite())
									currentMappings.add(ml);
							}
						}
						
						if(prevSize < currentMappings.size()) {
							minMissmatches = Integer.MAX_VALUE;
							for(MultiMappingLine ml : currentMappings) {
								if(ml.getMissmatches() < minMissmatches)
									minMissmatches = ml.getMissmatches();
							}
						}
					}
					
					Collections.sort(currentMappings);
					//in case we process a global context we have filter for reads which contain known splice signals or overlap known junctions
					currentMappings = getReadsOverlappingKnownSignals(currentMappings,minMissmatches);
					printRead(currentMappings,minMissmatches,this.maxMissmatchDifference,pw);
					
					minMissmatches = Integer.MAX_VALUE;
					if(missmatches < minMissmatches && (!isGlobalResolution || !isClippedLocation))
						minMissmatches = missmatches;
					
					foundCandidateSplit = isSplitLocation;
					
					prevReadId = currentReadId;
					currentMappings.clear();
					currentClippedMappings.clear();
					
					if(currentMappings.size() + currentClippedMappings.size() < mappingLinePool.size()) {
						tmpLine = mappingLinePool.get(currentMappings.size() + currentClippedMappings.size());
						tmpLine.update(currentLine,missmatches,-1,-1,isSplitLocation,hasKnownSpliceSignal,overlapsKnownJunction,overlapsCleavageSite);
					}
					else {
						tmpLine = new MultiMappingLine(currentLine,missmatches,-1,-1,isSplitLocation,hasKnownSpliceSignal,overlapsKnownJunction,overlapsCleavageSite);
						mappingLinePool.add(tmpLine);
					}
					
					if(!isGlobalResolution || !isClippedLocation) {
						currentMappings.add(tmpLine);
					}
					else {
						currentClippedMappings.add(tmpLine);
					}
					
					processedReads++;	
				}
			}
			
			//print the last read
			if(isGlobalResolution && !currentClippedMappings.isEmpty()) {
				if(currentMappings.isEmpty()) {
					currentMappings.addAll(currentClippedMappings);
					minMissmatches = Integer.MAX_VALUE;
					for(MultiMappingLine ml : currentMappings) {
						if(ml.getMissmatches() < minMissmatches)
							minMissmatches = ml.getMissmatches();
					}
				}
			}
			
			
			else if(!isGlobalResolution && !currentClippedMappings.isEmpty()) {
				//the situation in a default ContextMap run
				int prevSize = currentMappings.size();
				if(currentMappings.isEmpty()) {
					currentMappings.addAll(currentClippedMappings);
				}
				//clipped and regular alignments only appear in polyA tail recognition mode
				//if we found a candidate split alignment we will discard all clipped alignments.
				//otherwise we will add only clipped alignment overlapping with a cleavage site
				else if(!foundCandidateSplit) {
					for(MultiMappingLine ml : currentClippedMappings) {
						if(tmpLine.overlapsCleavageSite())
							currentMappings.add(ml);
					}
				}
				
				if(prevSize < currentMappings.size()) {
					minMissmatches = Integer.MAX_VALUE;
					for(MultiMappingLine ml : currentMappings) {
						if(ml.getMissmatches() < minMissmatches)
							minMissmatches = ml.getMissmatches();
					}
				}
			}
			
			Collections.sort(currentMappings);
			currentMappings = getReadsOverlappingKnownSignals(currentMappings,minMissmatches);
			printRead(currentMappings,minMissmatches,this.maxMissmatchDifference,pw);
			
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	private class SplitCandidateFilter extends Thread {
		private String inputPath;
		private String outputFilePath;
		
		private long startFilePointer;
		private long stopFilePointer;
		
		
		private int maxMismatchDifference;
		private int readLength;

		
		public SplitCandidateFilter(String inputPath, String outputFilePath, long startFilePointer, long stopFilePointer, int maxMismatchDifference, int readLength) {
			this.inputPath = inputPath;
			this.outputFilePath = outputFilePath;
			this.startFilePointer = startFilePointer;
			this.stopFilePointer = stopFilePointer;
			
			this.maxMismatchDifference = maxMissmatchDifference;
			this.readLength = readLength;
			
		}
		
		public void run() {
			try {
				
				BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputPath),"r",10000);
				braf.seek(this.startFilePointer);
				PrintWriter pw = new PrintWriter(new FileWriter(new File(this.outputFilePath)));
				String currentLine = braf.getNextLine();
				
				if(currentLine == null) {
					braf.close();
					return;
				}
				
				StringTokenizer st = new StringTokenizer(currentLine,"\t");
				String currentReadId = st.nextToken();
				String currentReadIdPrefix = currentReadId;
				HashSet<String> currentMscChrs = new HashSet<String>();
				boolean isMSC = false;
				if(currentReadId.contains("MSC")) {
					isMSC = true;
					currentReadIdPrefix = currentReadId.substring(0,currentReadId.indexOf("::"));
				}
				
				String prevReadId = currentReadId;
				String prevReadIdPrefix = currentReadIdPrefix;
				
				
				String mappingType = st.nextToken();
				//chr	start	end	strand 
				st.nextToken();
				int start = Integer.valueOf(st.nextToken());
				String endAsString = st.nextToken();
				int end;
				if(endAsString.equals("."))
				  end = start + readLength - 1;
				else
					end = Integer.valueOf(endAsString);
				
				String strand = st.nextToken();
				int mismatches = Integer.valueOf(st.nextToken());
				int minMismatches = Integer.MAX_VALUE;
				int mappingCount = 0;
				ArrayList<MultiMappingLine> currentSplitMappings = new ArrayList<MultiMappingLine>();
				ArrayList<MultiMappingLine> currentFullAndPartialMappings = new ArrayList<MultiMappingLine>();
				ArrayList<MultiMappingLine> currentMscMappings = new ArrayList<MultiMappingLine>();
				
				if(isMSC) {
					currentMscMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,false,'0',false,false));
					currentMscChrs.add(currentReadId.split("::")[2] + "_" + strand);
				}
				
				else if(mappingType.equals("S") && mismatches != Integer.MAX_VALUE) {
					currentSplitMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,true,'0',false,false));
				}
				else {
					currentFullAndPartialMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,false,'0',false,false));
					mappingCount++;
				}
				
				if(!isMSC && (mappingType.equals("S") || mappingType.equals("F")) && mismatches != Integer.MAX_VALUE) {
					minMismatches = mismatches;
				}
				
				while((currentLine = braf.getNextLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					currentReadId = st.nextToken();
					currentReadIdPrefix = currentReadId;
					isMSC = false;
					if(currentReadId.contains("MSC")) {
						isMSC = true;
						currentReadIdPrefix = currentReadId.substring(0,currentReadId.indexOf("::"));
					}
					
					mappingType = st.nextToken();
					//chr	start	end	strand 
					st.nextToken();
					start = Integer.valueOf(st.nextToken());
					endAsString = st.nextToken();
					if(endAsString.equals("."))
					  end = start + readLength - 1;
					else
						end = Integer.valueOf(endAsString);
					
					strand = st.nextToken();
					mismatches = Integer.valueOf(st.nextToken());
					
					if(prevReadId.equals(currentReadId) || prevReadIdPrefix.equals(currentReadIdPrefix)) {
						if(isMSC) {
							currentMscChrs.add(currentReadId.split("::")[2] + "_" + strand);
							currentMscMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,false,'0',false,false));
						}
						
						else if(mappingType.equals("S") && mismatches != Integer.MAX_VALUE) {
							currentSplitMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,true,'0',false,false));
						}
						else {
							currentFullAndPartialMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,false,'0',false,false));
							mappingCount++;
						}
						
						if(!isMSC && (mappingType.equals("S") || mappingType.equals("F")) && mismatches != Integer.MAX_VALUE && mismatches < minMismatches)
							minMismatches = mismatches;
					}
					
					else {
						
						mappingCount += currentMscChrs.size() + getNumberOfPrintableSplitReads(currentSplitMappings,minMismatches,maxMissmatchDifference);
						
						//print all full and partial reads
						for(MultiMappingLine l : currentFullAndPartialMappings) {
							pw.println(l.getLine() + "\t" + mappingCount + "\t0");
						}
						
						//print all msc reads
						for(MultiMappingLine l : currentMscMappings) {
							pw.println(l.getLine() + "\t" + mappingCount + "\t0");
						}
						
						//print splits
						printRead(currentSplitMappings,minMismatches,maxMissmatchDifference,mappingCount,pw);
						
						
						prevReadId = currentReadId;
						prevReadIdPrefix = currentReadIdPrefix;
						
						currentMscMappings.clear();
						currentMscChrs.clear();
						currentFullAndPartialMappings.clear();
						currentSplitMappings.clear();
						mappingCount = 0;
						
						if(isMSC) {
							currentMscMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,false,'0',false,false));
							currentMscChrs.add(currentReadId.split("::")[2] + "_" + strand);
						}
						
						else if(mappingType.equals("S") && mismatches != Integer.MAX_VALUE) {
							currentSplitMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,true,'0',false,false));
						}
						
						else {
							currentFullAndPartialMappings.add(new MultiMappingLine(currentLine,mismatches,start,end,false,'0',false,false));
							mappingCount++;
						}
						
						minMismatches = Integer.MAX_VALUE;
						if(!isMSC && (mappingType.equals("S") || mappingType.equals("F")) && mismatches != Integer.MAX_VALUE && mismatches < minMismatches)
							minMismatches = mismatches;
					}
					
					
					if(braf.getFilePointer() == this.stopFilePointer)
						break;
					
				}
				
				mappingCount += currentMscChrs.size();
				//print all full and partial reads
				for(MultiMappingLine l : currentFullAndPartialMappings) {
					pw.println(l.getLine() + "\t" + mappingCount + "\t0");
				}
				
				//print all msc reads
				for(MultiMappingLine l : currentMscMappings) {
					pw.println(l.getLine() + "\t" + mappingCount + "\t0");
				}
			
				//print the last split reads
				printRead(currentSplitMappings,minMismatches,maxMissmatchDifference,mappingCount,pw);
				
				braf.close();
				pw.close();
				
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
		
		private int getNumberOfPrintableSplitReads(ArrayList<MultiMappingLine> mappings, int minMissmatches,int maxMissmatchDifference) {
			int counter = 0;
			for(int i = 0; i < mappings.size(); i++) {
				if(mappings.get(i).getMissmatches() - minMissmatches <= maxMissmatchDifference) {
					counter++;
				}
			}
			return counter;
		}
		
		
		private void printRead(ArrayList<MultiMappingLine> mappings, int minMissmatches,int maxMissmatchDifference, int mappingCount, PrintWriter pw) throws Exception {
			for(int i = 0; i < mappings.size(); i++) {
				if(mappings.get(i).getMissmatches() - minMissmatches <= maxMissmatchDifference) {
					pw.println(mappings.get(i).getLine() + "\t" + mappingCount + "\t0");
				}
			}
		}
	}
	
	
	public void filterSplitCandidates(String inputPath, String outputDirPath,String outputFilePath, int maxMismatchDifference, int readLength, int threads, boolean pairedEnd) {
		try {
			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputPath),"r",1024);
			long fileSize = braf.getChannel().size();
			long incRate = fileSize/threads;
			long prevPosition = 0;
			long currentPosition = 0;
			String tmpOutputFilePath;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readId;
			String prevReadId;
			int chunkIndex = 0;
			
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				braf.seek(currentPosition);
				braf.getNextLine();
				currentLine = braf.getNextLine();
				prevReadId = null;
				currentPosition = braf.getFilePointer();
				while(currentLine != null) {
					splittedLine = tabPattern.split(currentLine);
					readId = splittedLine[0];
					if(pairedEnd) {
						readId = splittedLine[0].substring(0,splittedLine[0].lastIndexOf("/"));
						if(splittedLine[0].contains("MSC"))
							readId = splittedLine[0].substring(0,splittedLine[0].indexOf("::") - 2);
					}
					if(prevReadId == null) {
						prevReadId = readId;
					}
					else if(!readId.equals(prevReadId)) {
						break;
					}
					
					currentPosition = braf.getFilePointer();
					currentLine = braf.getNextLine();
				}
				
				tmpOutputFilePath = outputDirPath + "/filtered_splits_" + chunkIndex + ".rmap";
				futures.add(executor.submit(new SplitCandidateFilter(inputPath, tmpOutputFilePath, prevPosition, currentPosition, maxMismatchDifference, readLength)));
				
				prevPosition = currentPosition;
				chunkIndex++;
			}
			braf.close();
			
			executor.shutdown();
			for(Future f : futures) {
				f.get();
			}
			
			concatenateFiles(outputDirPath,outputFilePath);
			File[] tmpFiles = new File(outputDirPath).listFiles();
			for(File f : tmpFiles)
				f.delete();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	private void concatenateFiles(String inputDirPath, String outputFilePath) throws Exception {
		File outputFile = new File(outputFilePath);
		if(outputFile.isFile())
			outputFile.delete();
		
		FileOutputStream fos = new FileOutputStream(outputFilePath,true);
		FileChannel writeChannel = fos.getChannel();
		RandomAccessFile rf;
		FileChannel readChannel;
		File[] tmpFiles = new File(inputDirPath).listFiles();
		long currentChannelSize;
		long transferedBytes;
		for(File tmpFile : tmpFiles) {
			rf = new RandomAccessFile(tmpFile,"r");
			readChannel = rf.getChannel();
			currentChannelSize = readChannel.size();
			transferedBytes = 0;
			while(transferedBytes < currentChannelSize) {
				transferedBytes += readChannel.transferTo(transferedBytes, readChannel.size(), writeChannel);
			}
			rf.close();
		}
		fos.close();
	}
	
	
	//Expects a multi mapping file which is sorted by read id.
	//
/*	public void filterSplitCandidates(String inputPath, String outputPath,int maxMissmatchDifference, int readLength) {
		try {
			
			BufferedReader br = new BufferedReader(new FileReader(new File(inputPath)));
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(outputPath));			
			String currentLine = br.readLine();
			
			if(currentLine == null) {
				br.close();
				pw.close();
				return;
			}
			
			StringTokenizer st = new StringTokenizer(currentLine,"\t");
			String currentReadId = st.nextToken();
			String prevReadId = currentReadId;
			String mappingType = st.nextToken();
			//chr	start	end	strand 
			st.nextToken();
			int start = Integer.valueOf(st.nextToken());
			String endAsString = st.nextToken();
			int end;
			if(endAsString.equals("."))
			  end = start + readLength - 1;
			else
				end = Integer.valueOf(endAsString);
			
			st.nextToken();
			int missmatches = Integer.valueOf(st.nextToken());
			int minMissmatches = Integer.MAX_VALUE;
			ArrayList<MultiMappingLine> currentMappings = new ArrayList<MultiMappingLine>();
			if(mappingType.equals("S") && missmatches != Integer.MAX_VALUE) {
				currentMappings.add(new MultiMappingLine(currentLine,missmatches,start,end,true,false,false));
			}
			else {
				pw.write(currentLine);
				pw.newLine();
			}
			
			if((mappingType.equals("S") || mappingType.equals("F")) && missmatches != Integer.MAX_VALUE) {
				minMissmatches = missmatches;
			}
			
			while(br.ready()) {
				currentLine = br.readLine();
				st = new StringTokenizer(currentLine,"\t");
				currentReadId = st.nextToken();
				mappingType = st.nextToken();
				//chr	start	end	strand 
				st.nextToken();
				start = Integer.valueOf(st.nextToken());
				endAsString = st.nextToken();
				if(endAsString.equals("."))
				  end = start + readLength - 1;
				else
					end = Integer.valueOf(endAsString);
				st.nextToken();
				missmatches = Integer.valueOf(st.nextToken());
				
				if(prevReadId.equals(currentReadId)) {
					if(mappingType.equals("S") && missmatches != Integer.MAX_VALUE) {
						currentMappings.add(new MultiMappingLine(currentLine,missmatches,start,end,true,false,false));
					}
					else {
						pw.write(currentLine);
						pw.newLine();
					}
					
					if((mappingType.equals("S") || mappingType.equals("F")) && missmatches != Integer.MAX_VALUE && missmatches < minMissmatches)
						minMissmatches = missmatches;
				}
				
				else {
					Collections.sort(currentMappings);
					//print splits
					printRead(currentMappings,minMissmatches,maxMissmatchDifference,pw);
					
					prevReadId = currentReadId;
					currentMappings.clear();
					if(mappingType.equals("S") && missmatches != Integer.MAX_VALUE) {
						currentMappings.add(new MultiMappingLine(currentLine,missmatches,start,end,true,false,false));
					}
					
					else {
						pw.write(currentLine);
						pw.newLine();
					}
					
					minMissmatches = Integer.MAX_VALUE;
					if((mappingType.equals("S") || mappingType.equals("F")) && missmatches != Integer.MAX_VALUE && missmatches < minMissmatches)
						minMissmatches = missmatches;
				}
			}
			//print the last split reads
			printRead(currentMappings,minMissmatches,maxMissmatchDifference,pw);
			
			br.close();
			pw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
*/	
	
	
	private class PairedEndFilter extends Thread {
		
		private String inputFilePath;
		private String outputFilePath;
		private int readLength;
		private int maxContextSize;
		
		private long startPointer;
		private long stopPointer;
		
		
		public PairedEndFilter(String inputFilePath, long startPointer, long stopPointer, String outputFilePath, int readLength, int maxContextSize) {
			this.inputFilePath = inputFilePath;
			this.startPointer = startPointer;
			this.stopPointer = stopPointer;
			this.outputFilePath = outputFilePath;
			this.readLength = readLength;
			this.maxContextSize = maxContextSize;
		}
		
		public void run() {
			try {
				BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputFilePath), "r",10000);
				UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(outputFilePath));
				
				braf.seek(this.startPointer);			
				String currentLine = braf.getNextLine();
				if(currentLine == null) {
					braf.close();
					pw.close();
					return;
				}
				
				
				StringTokenizer st = new StringTokenizer(currentLine,"\t");
				String readId = st.nextToken();
				String readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
				boolean isMSC = false;
				if(readId.contains("MSC")) {
					readIdPrefix = readId.substring(0,readId.indexOf("::") - 2);
					isMSC = true;
					
				}
				
				
				int mateNumber = Integer.valueOf(readId.substring(readId.length() - 1));
				String prevReadId = readId;
				String prevReadIdPrefix = readIdPrefix;
				char mappingType = st.nextToken().charAt(0);
				String chr = st.nextToken(); 
				int start = Integer.valueOf(st.nextToken());
				String endAsString = st.nextToken();
				char strand = st.nextToken().charAt(0);
				int mismatches = Integer.valueOf(st.nextToken());
				if(st.hasMoreTokens()) {
					readLength = Integer.valueOf(st.nextToken());
				}
				
				
				int end;
				if(endAsString.equals("."))
				  end = start + readLength - 1;
				else
					end = Integer.valueOf(endAsString);
				
				Read currentRead = new Read(readId);
				ReadPair<Read,Read> currentPair = new ReadPair<Read,Read>(readIdPrefix);
				currentPair.setFirst(null);
				currentPair.setSecond(null);
				
				Pair<ArrayList<String>,ArrayList<String>> currentMappingLines = new Pair<ArrayList<String>,ArrayList<String>>();
				currentMappingLines.setFirst(new ArrayList<String>());
				currentMappingLines.setSecond(new ArrayList<String>());
				ArrayList<String> currentMappingList;
				
				if(mateNumber == 1) {
					currentPair.setFirst(currentRead);
					currentMappingList = currentMappingLines.getFirst();
				
				}
				else {
					currentPair.setSecond(currentRead);
					currentMappingList = currentMappingLines.getSecond();
				
				}
				
				ReadLocation newLocation = new ReadLocation(chr, strand, mappingType, start, end, -1, -1, mismatches);
				newLocation.setReadId(readId);
				newLocation.setMsc(isMSC);
				currentRead.addLocation(newLocation);
				currentMappingList.add(currentLine);

				PairedEndPair validPairs;
				String[] splittedLine;
				Pattern tabPattern = Pattern.compile("\t");
				
				boolean mate1hasPartialAlignmentsOnly;
				boolean mate2hasPartialAlignmentsOnly;
				
				while((currentLine = braf.getNextLine()) != null) {
					st = new StringTokenizer(currentLine,"\t");
					readId = st.nextToken();
					readIdPrefix = readId.substring(0,readId.lastIndexOf("/"));
					isMSC = false;
					if(readId.contains("MSC")) {
						readIdPrefix = readId.substring(0,readId.indexOf("::") - 2);
						isMSC = true;
					}
					
					mateNumber = Integer.valueOf(readId.substring(readId.length() - 1));
					mappingType = st.nextToken().charAt(0);
					chr = st.nextToken(); 
					start = Integer.valueOf(st.nextToken());
					endAsString = st.nextToken();
					strand = st.nextToken().charAt(0);
					mismatches = Integer.valueOf(st.nextToken());
					if(st.hasMoreTokens()) {
						readLength = Integer.valueOf(st.nextToken());
					}
					
					if(endAsString.equals("."))
						  end = start + readLength - 1;
						else
							end = Integer.valueOf(endAsString);
					
				
					
					
					if(readId.equals(prevReadId)) {
							newLocation = new ReadLocation(chr, strand, mappingType, start, end, -1,-1, mismatches);
							newLocation.setReadId(readId);
							newLocation.setMsc(isMSC);
							currentRead.addLocation(newLocation);
							currentMappingList.add(currentLine);
					}
					
					//second read from the actual fragment
					else if(readIdPrefix.equals(prevReadIdPrefix)) {
						
						if(mateNumber == 1 && currentPair.getFirst() != null)
								currentRead = currentPair.getFirst();
						
						else if(mateNumber == 2 && currentPair.getSecond() != null)
								currentRead = currentPair.getSecond();
						
						else
							currentRead = new Read(readId);
							
						
						if(mateNumber == 1) {
							currentPair.setFirst(currentRead);
							currentMappingList = currentMappingLines.getFirst();
							
						}
						else {
							currentPair.setSecond(currentRead);
							currentMappingList = currentMappingLines.getSecond();
							
						}
						
						newLocation = new ReadLocation(chr, strand, mappingType, start, end, -1,-1, mismatches);
						newLocation.setReadId(readId);
						newLocation.setMsc(isMSC);
						currentRead.addLocation(newLocation);
						currentMappingList.add(currentLine);
						prevReadId = readId;
					}
					
					//here we have a new fragment
					else {
						
						validPairs = getReadsBuildingValidPairs(currentPair,maxContextSize);
						
						//no valid pair found
						if(validPairs.getFirst().isEmpty()) {
							
							//only print out alignments iff not both mates have alignments
							//NOTE: We change this now and wait until the resolution of the local contexts if one of the mates gets discarded completely.
							//This can only happen, if it has partial alignments only
							mate1hasPartialAlignmentsOnly = true;
							mate2hasPartialAlignmentsOnly = true;
							if(currentPair.getFirst() != null && currentPair.getSecond() != null) {
								for(ReadLocation l : currentPair.getFirst().getLocations()) {
									if(l.getMappingType() == 'F' || l.getMappingType() == 'S') {
										mate1hasPartialAlignmentsOnly = false;
										break;
									}
								}
								for(ReadLocation l : currentPair.getSecond().getLocations()) {
									if(l.getMappingType() == 'F' || l.getMappingType() == 'S') {
										mate2hasPartialAlignmentsOnly = false;
										break;
									}
								}
							}
							
							if(currentPair.getFirst() == null || currentPair.getSecond() == null || mate1hasPartialAlignmentsOnly || mate2hasPartialAlignmentsOnly) {
								for(String line : currentMappingLines.getFirst()) {
									pw.write(line);
									pw.newLine();
								}
								
								
								for(String line : currentMappingLines.getSecond()) {
									pw.write(line);
									pw.newLine();
								}
							}
						}
						
						//at least one valid pair found -> only print out the reads part of a valid pair
						else {
							for(int i : validPairs.getFirst()) {
								splittedLine = tabPattern.split(currentMappingLines.getFirst().get(i));
								splittedLine[splittedLine.length - 2] = Integer.toString(validPairs.getMappingCountA());
								splittedLine[splittedLine.length - 1] = Integer.toString(validPairs.getValidPairCount());
								pw.write(splittedLine[0]);
								for(int j = 1; j < splittedLine.length; j++) {
									pw.write("\t" + splittedLine[j]);
								}
								pw.newLine();
							}
							
							for(int i : validPairs.getSecond()) {
								splittedLine = tabPattern.split(currentMappingLines.getSecond().get(i));
								splittedLine[splittedLine.length - 2] = Integer.toString(validPairs.getMappingCountB());
								splittedLine[splittedLine.length - 1] = Integer.toString(validPairs.getValidPairCount());
								
								pw.write(splittedLine[0]);
								for(int j = 1; j < splittedLine.length; j++) {
									pw.write("\t" + splittedLine[j]);
								}
								pw.newLine();
							}
						}
							
						
						prevReadId = readId;
						prevReadIdPrefix = readIdPrefix;
						
						currentRead = new Read(readId);
						currentPair = new ReadPair<Read,Read>(readIdPrefix);
						currentPair.setFirst(null);
						currentPair.setSecond(null);
						
						currentMappingLines.getFirst().clear();
						currentMappingLines.getSecond().clear();
						
						if(mateNumber == 1) {
							currentPair.setFirst(currentRead);
							currentMappingList = currentMappingLines.getFirst();
						}
						else {
							currentPair.setSecond(currentRead);
							currentMappingList = currentMappingLines.getSecond();
						}

						newLocation = new ReadLocation(chr, strand, mappingType, start, end, -1, -1, mismatches);
						newLocation.setReadId(readId);
						newLocation.setMsc(isMSC);
						currentRead.addLocation(newLocation);
						currentMappingList.add(currentLine);
						
					}
				
					if(braf.getFilePointer() == this.stopPointer)
						break;
				}
				
				
				//process last lines
				validPairs = getReadsBuildingValidPairs(currentPair,maxContextSize);
				
				//no valid pair found
				if(validPairs.getFirst().isEmpty()) {
					
					mate1hasPartialAlignmentsOnly = true;
					mate2hasPartialAlignmentsOnly = true;
					if(currentPair.getFirst() != null && currentPair.getSecond() != null) {
						for(ReadLocation l : currentPair.getFirst().getLocations()) {
							if(l.getMappingType() == 'F' || l.getMappingType() == 'S') {
								mate1hasPartialAlignmentsOnly = false;
								break;
							}
						}
						for(ReadLocation l : currentPair.getSecond().getLocations()) {
							if(l.getMappingType() == 'F' || l.getMappingType() == 'S') {
								mate2hasPartialAlignmentsOnly = false;
								break;
							}
						}
					}
					
					if(currentPair.getFirst() == null || currentPair.getSecond() == null || mate1hasPartialAlignmentsOnly || mate2hasPartialAlignmentsOnly) {
						for(String line : currentMappingLines.getFirst()) {
							pw.write(line);
							pw.newLine();
						}
						
						
						for(String line : currentMappingLines.getSecond()) {
							pw.write(line);
							pw.newLine();
						}
					}
				}
				
				//at least one valid pair found -> only print out the reads part of a valid pair
				else {
					for(int i : validPairs.getFirst()) {
						splittedLine = tabPattern.split(currentMappingLines.getFirst().get(i));
						splittedLine[splittedLine.length - 2] = Integer.toString(validPairs.getMappingCountA());
						splittedLine[splittedLine.length - 1] = Integer.toString(validPairs.getValidPairCount());
						
						
						pw.write(splittedLine[0]);
						for(int j = 1; j < splittedLine.length; j++) {
							pw.write("\t" + splittedLine[j]);
						}
						pw.newLine();
					}
					
					for(int i : validPairs.getSecond()) {
						splittedLine = tabPattern.split(currentMappingLines.getSecond().get(i));
						splittedLine[splittedLine.length - 2] = Integer.toString(validPairs.getMappingCountB());
						splittedLine[splittedLine.length - 1] = Integer.toString(validPairs.getValidPairCount());
						
						pw.write(splittedLine[0]);
						for(int j = 1; j < splittedLine.length; j++) {
							pw.write("\t" + splittedLine[j]);
						}
						pw.newLine();
					}
				}
				
				
				braf.close();
				pw.close();
			}
			
			catch(Exception e) {
				e.printStackTrace();
			}
		}
		
	}
	
	public void filterPairedEndCandidates(String inputFilePath, String outputDirPath, String outputFilePath, int readLength, int maxContextSize, int threads) {
		try {
			BufferedRandomAccessFile braf = new BufferedRandomAccessFile(new File(inputFilePath),"r",1024);
			long fileSize = braf.getChannel().size();
			long incRate = fileSize/threads;
			long prevPosition = 0;
			long currentPosition = 0;
			String tmpOutputFilePath;
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			String readIdPrefix;
			String prevReadIdPrefix;
			int chunkIndex = 0;
			
			ExecutorService executor = Executors.newFixedThreadPool(threads);
			ArrayList<Future> futures = new ArrayList<Future>();
			
			while(currentPosition < fileSize) {
				currentPosition += incRate;
				braf.seek(currentPosition);
				braf.getNextLine();
				currentLine = braf.getNextLine();
				
				prevReadIdPrefix = null;
				currentPosition = braf.getFilePointer();
				while(currentLine != null) {
					splittedLine = tabPattern.split(currentLine);
					readIdPrefix = splittedLine[0].substring(0,splittedLine[0].lastIndexOf("/"));
					if(splittedLine[0].contains("MSC"))
						readIdPrefix = splittedLine[0].substring(0,splittedLine[0].indexOf("::") - 2);
					if(prevReadIdPrefix == null) {
						prevReadIdPrefix = readIdPrefix;
					}
					else if(!readIdPrefix.equals(prevReadIdPrefix)) {
						break;
					}
					
					currentPosition = braf.getFilePointer();
					currentLine = braf.getNextLine();
				}
				
				tmpOutputFilePath = outputDirPath + "/filtered_mates_" + chunkIndex + ".rmap";
				futures.add(executor.submit(new PairedEndFilter(inputFilePath, prevPosition, currentPosition, tmpOutputFilePath, readLength, maxContextSize)));
				
				prevPosition = currentPosition;
				chunkIndex++;
			}
			braf.close();
			
			executor.shutdown();
			for(Future f : futures) {
				f.get();
			}
			
			concatenateFiles(outputDirPath,outputFilePath);
			File[] tmpFiles = new File(outputDirPath).listFiles();
			for(File f : tmpFiles)
				f.delete();
		
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * currently this function only uses the following constraints:
	 * 
	 * - mates have to be aligned to opposite strands (first: + -> second: - or first: - -> second: +)
	 * - the distance between two mates is not allowed to be larger than the maximum context size (see variable maxContextSize).
	 *
	 * 
	 * TODO: calculate a more accurate fragment length distribution from fully aligned reads to further restrict the distance criterion.
	 * TODO: Think about a maximum allowed context-score difference between two mates. It makes no sense to considers cases where one mate
	 * has almost no support by other reads whereas the other one is seems to originate from a highly expressed region.
	 * 
	 * 
	 *
	 * 
	 **/
	
	private PairedEndPair getReadsBuildingValidPairs(ReadPair<Read,Read> readPair, int maxContextSize) {
		try {
			
			int pairCount = 0;
			HashSet<Integer> mappingCountA = new HashSet<Integer>();
			HashSet<Integer> mappingCountB = new HashSet<Integer>();
			PairedEndPair validPairs = new PairedEndPair();
			validPairs.setFirst(new HashSet<Integer>());
			validPairs.setSecond(new HashSet<Integer>());
			
			if(readPair.getFirst() == null || readPair.getSecond() == null) {
				validPairs.setValidPairCount(0);
				return validPairs;
			}
			
			ArrayList<ReadLocation> locationsA;
			ArrayList<ReadLocation> locationsB;
			locationsA = readPair.getFirst().getLocations();
			int currentEnd;

			HashSet<String> mscA = new HashSet<String>();
			HashSet<String> mscB = new HashSet<String>();
			
			locationsB = readPair.getSecond().getLocations();
			for(int i = 0; i < locationsA.size(); i++) {
				mscB.clear();
				for(int j = 0; j < locationsB.size(); j++) {
					
					//chromosome criterion not fulfilled
					if(!locationsA.get(i).getChr().equals(locationsB.get(j).getChr()))
						continue;
						
					//strand criterion not fulfilled
					if(locationsA.get(i).getStrand() == locationsB.get(j).getStrand())
						continue;
						
					
					//A starts upstream of B
					if(locationsA.get(i).getStartA() <= locationsB.get(j).getStartA()) {
						//distance criterion fulfilled
						currentEnd = locationsA.get(i).getEndA();
						
						
						if(Math.abs(currentEnd - locationsB.get(j).getStartA()) <= maxContextSize) {
							validPairs.getFirst().add(i);
							validPairs.getSecond().add(j);
							
							if((!locationsA.get(i).isMSC() || !mscA.contains(locationsA.get(i).getChr() + "_" + locationsA.get(i).getStrand())) &&
							    (!locationsB.get(j).isMSC() || !mscB.contains(locationsB.get(j).getChr() + "_" + locationsB.get(j).getStrand()))) {	
								pairCount++;
								mappingCountA.add(i);
								mappingCountB.add(j);
								
								if(locationsB.get(j).isMSC())
									mscB.add(locationsB.get(j).getChr() + "_" + locationsB.get(j).getStrand());
							}
						}
					}
					
					//B starts upstream of A
					else if(locationsA.get(i).getStartA() > locationsB.get(j).getStartA()) {
						
						//distance criterion fulfilled
						currentEnd = locationsB.get(j).getEndA();
						
						if(Math.abs(currentEnd - locationsA.get(i).getStartA()) <= maxContextSize) {
							validPairs.getFirst().add(i);
							validPairs.getSecond().add(j);
							
							if((!locationsA.get(i).isMSC() || !mscA.contains(locationsA.get(i).getChr() + "_" + locationsA.get(i).getStrand())) &&
								    (!locationsB.get(j).isMSC() || !mscB.contains(locationsB.get(j).getChr() + "_" + locationsB.get(j).getStrand()))) {	
									pairCount++;
									mappingCountA.add(i);
									mappingCountB.add(j);
									
									if(locationsB.get(j).isMSC())
										mscB.add(locationsB.get(j).getChr() + "_" + locationsB.get(j).getStrand());
							}
						}
					}
					
					
				}
				
				if(locationsA.get(i).isMSC() && validPairs.getFirst().contains(i)) {
					mscA.add(locationsA.get(i).getChr() + "_" + locationsA.get(i).getStrand());
				}
				
			}
			
			validPairs.setValidPairCount(pairCount);
			validPairs.setMappingCountA(mappingCountA.size());
			validPairs.setMappingCountB(mappingCountB.size());
		
			return validPairs;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	
	private class PairedEndPair extends Pair<HashSet<Integer>,HashSet<Integer>> {
		protected int validPairCount;
		protected int mappingCountA;
		protected int mappingCountB;
		
		public PairedEndPair() {
			super();
			this.validPairCount = 0;
			this.mappingCountA = 0;
			this.mappingCountB = 0;
		}
		
		public int getValidPairCount() {
			return this.validPairCount;
		}
		
		public void setValidPairCount(int count) {
			this.validPairCount = count;
		}
		
		public int getMappingCountA() {
			return this.mappingCountA;
		}
		
		public void setMappingCountA(int count) {
			this.mappingCountA = count;
		}
		
		public int getMappingCountB() {
			return this.mappingCountB;
		}
		
		public void setMappingCountB(int count) {
			this.mappingCountB = count;
		}
	}
	
	
	
	/**
	 * currently this function only uses the following constraints:
	 * 
	 * - mates have to be aligned to opposite strands (first: + -> second: - or first: - -> second: +)
	 * - the distance between two mates is not allowed to be larger than the maximum context size (see variable maxContextSize).
	 *
	 * 
	 * TODO: calculate a more accurate fragment length distribution from fully aligned reads to further restrict the distance criterion.
	 * TODO: Think about a maximum allowed context-score difference between two mates. It makes no sense to considers cases where one mate
	 * has almost no support by other reads whereas the other one is seems to originate from a highly expressed region.
	 * 
	 * 
	 *
	 * 
	 **/
	
	private void getValidReadPairs(ReadPair<Read,Read> readPair,int maxContextSize) {
		try {
			ArrayList<ReadLocation> locationsA;
			ArrayList<ReadLocation> locationsB;
			locationsA = readPair.getFirst().getLocations();
		
			int currentStart;
			int currentEnd;
			int currentDistance;
			int secondLocationEnd;
			int coordinatesIndexBeforeClippingA;
			int coordinatesIndexBeforeClippingB;
			if(readPair.getSecond() != null) {
				locationsB = readPair.getSecond().getLocations();
				for(int i = 0; i < locationsA.size(); i++) {
					
					coordinatesIndexBeforeClippingA = locationsA.get(i).getCoordinates().size() - 1;
					for(int k = 0; k < locationsA.get(i).getCoordinates().size(); k++) {
						Pair<Integer,Integer> segment = locationsA.get(i).getCoordinates().get(k);
						
						if(segment.getFirst() <= 0) {
							coordinatesIndexBeforeClippingA = k - 1;
							break;
						}
					}
					
					for(int j = 0; j < locationsB.size(); j++) {
						
						if(!locationsA.get(i).getChr().equals(locationsB.get(j).getChr()))
							continue;
						
						coordinatesIndexBeforeClippingB = locationsB.get(j).getCoordinates().size() - 1;
						for(int k = 0; k < locationsB.get(j).getCoordinates().size(); k++) {
							Pair<Integer,Integer> segment = locationsB.get(j).getCoordinates().get(k);
							
							if(segment.getFirst() <= 0) {
								coordinatesIndexBeforeClippingB = k - 1;
								break;
							}
						}
						
						//strand criterion fulfilled
						if(locationsA.get(i).getStrand() != locationsB.get(j).getStrand()) {
							
							//A starts upstream of B
							if(locationsA.get(i).getCoordinates().get(0).getFirst() <= locationsB.get(j).getCoordinates().get(0).getFirst()) {
								//distance criterion fulfilled
								currentStart = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getFirst();
								currentEnd = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getSecond();
/*								if(locationsA.get(i).getCoordinates().size() > 1 && locationsA.get(i).getCoordinates().get(1).getFirst() > 0) {
									currentStart = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getFirst();
									currentEnd = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size()-1).getSecond();
								}
*/								
								currentDistance = locationsB.get(j).getCoordinates().get(0).getFirst() - currentEnd;
								
								//distance criterion fulfilled
								if(currentDistance >= 0 && currentDistance <= maxContextSize) {
									readPair.addValidPair(new Pair<Integer,Integer>(i,j));
								}
								else if(currentDistance < 0) {
									secondLocationEnd = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getSecond();
/*									if(locationsB.get(j).getCoordinates().size() > 1 && locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size() - 1).getSecond() > 0)
										secondLocationEnd = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size() - 1).getSecond();
*/									
									if(currentEnd < secondLocationEnd && secondLocationEnd - currentEnd <= maxContextSize) {
										readPair.addValidPair(new Pair<Integer,Integer>(i,j));
									}
								}
							}
							
							//B starts upstream of A
							else if(locationsA.get(i).getCoordinates().get(0).getFirst() > locationsB.get(j).getCoordinates().get(0).getFirst()) {
								//distance criterion fulfilled
								currentStart = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getFirst();
								currentEnd = locationsB.get(j).getCoordinates().get(coordinatesIndexBeforeClippingB).getSecond();
/*								if(locationsB.get(j).getCoordinates().size() > 1 && locationsB.get(j).getCoordinates().get(1).getFirst() > 0) {
									currentStart = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getFirst();
									currentEnd = locationsB.get(j).getCoordinates().get(locationsB.get(j).getCoordinates().size()-1).getSecond();
								}
*/								
								currentDistance = locationsA.get(i).getCoordinates().get(0).getFirst() - currentEnd;
								if(currentDistance >= 0 && currentDistance <= maxContextSize) {
									readPair.addValidPair(new Pair<Integer,Integer>(i,j));
								}
								else if(currentDistance < 0) {
									secondLocationEnd = locationsA.get(i).getCoordinates().get(coordinatesIndexBeforeClippingA).getSecond();
/*									if(locationsA.get(i).getCoordinates().size() > 1 && locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size() - 1).getSecond() > 0)
										secondLocationEnd = locationsA.get(i).getCoordinates().get(locationsA.get(i).getCoordinates().size() - 1).getSecond();
*/									
									if(currentEnd < secondLocationEnd && secondLocationEnd - currentEnd <= maxContextSize) {
										readPair.addValidPair(new Pair<Integer,Integer>(i,j));
									}
								}
							}
						}
					}
				}
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	private ArrayList<MultiMappingLine> getReadsOverlappingKnownSignals(ArrayList<MultiMappingLine> mappings, int minMissmatches) {
		ArrayList<MultiMappingLine> filteredMappings = new ArrayList<MultiMappingLine>();
		boolean foundMappingWithKnownSpliceSignal = false;
		boolean foundMappingOverlappingKnownJunction = false;
		for(int i  = 0;i< mappings.size(); i++) {
			if(mappings.get(i).getMissmatches() - minMissmatches <= this.maxMissmatchDifference) {
				if(mappings.get(i).getStrandOfSpliceSignal() != '0') foundMappingWithKnownSpliceSignal = true;
				if(mappings.get(i).overlapsKnownJunction()) foundMappingOverlappingKnownJunction = true;
			}
			else
				break;
		}
		
		for(int i  = 0;i< mappings.size(); i++) {
			if(mappings.get(i).getMissmatches() - minMissmatches <= this.maxMissmatchDifference) {
				if(!mappings.get(i).isSplitLocation()) {
					filteredMappings.add(mappings.get(i));
					continue;
				}
				
				if(!foundMappingOverlappingKnownJunction && (!foundMappingWithKnownSpliceSignal || !this.preferExtensionsWithKnownSpliceSignal)) {
					filteredMappings.add(mappings.get(i));
					continue;
				}
				
				if(foundMappingOverlappingKnownJunction && mappings.get(i).overlapsKnownJunction()) {
					filteredMappings.add(mappings.get(i));
					continue;
				}
				
				if(!foundMappingOverlappingKnownJunction && foundMappingWithKnownSpliceSignal && mappings.get(i).getStrandOfSpliceSignal() != '0') {
					filteredMappings.add(mappings.get(i));
					continue;
				}
					
			}
			else
				break;
		}
		
		return filteredMappings;
	}
	
	private void printRead(ArrayList<MultiMappingLine> mappings, int minMissmatches,int maxMissmatchDifference, UnsynchronizedBufferedWriter pw) throws Exception {
			for(int i = 0; i < mappings.size(); i++) {
				if(mappings.get(i).getMissmatches() - minMissmatches <= maxMissmatchDifference) {
					pw.write(mappings.get(i).getLine());
					pw.newLine();
				}
				else
					break;
			}
	}
	
/*	private void printPartialReads(ArrayList<MultiMappingLine> mappings, IntervalTree<PartialReadAlignment> partialReads, int minMissmatches,int maxMissmatchDifference, UnsynchronizedBufferedWriter pw) throws Exception {
		Collection<PartialReadAlignment> currentPartialReads;
		for(int i = 0; i < mappings.size(); i++) {
			if(mappings.get(i).getMissmatches() - minMissmatches <= maxMissmatchDifference) {
				currentPartialReads = partialReads.getIntervalsSpanning(mappings.get(i).getStart(), new ArrayList<PartialReadAlignment>());
				for(PartialReadAlignment tmpRead : currentPartialReads) {
					if(tmpRead.getStart() == mappings.get(i).getStart())
						partialReads.remove(tmpRead);
				}
				
				currentPartialReads = partialReads.getIntervalsSpanning(mappings.get(i).getEnd(), new ArrayList<PartialReadAlignment>());
				for(PartialReadAlignment tmpRead : currentPartialReads) {
					if(tmpRead.getStop() == mappings.get(i).getEnd())
						partialReads.remove(tmpRead);
				}
			}
			else
				break;
		}
		Iterator<PartialReadAlignment> readIt = partialReads.iterator();
		while(readIt.hasNext()) {
			pw.write(readIt.next().getLine());
			pw.newLine();
		}
	}
	
	private class PartialReadAlignment implements Interval {

		private int start;
		private int stop;
		private String alignmentLine;
		
		public PartialReadAlignment(int start, int stop, String alignmentLine) {
			this.start = start;
			this.stop = stop;
			this.alignmentLine = alignmentLine;
		}
		
		public String getLine() {
			return this.alignmentLine;
		}
		
		@Override
		public int getStart() {
			// TODO Auto-generated method stub
			return this.start;
		}

		@Override
		public int getStop() {
			// TODO Auto-generated method stub
			return this.stop;
		}
		
	}*/
	
	
	private class MultiMappingLine implements Comparable {

		private String mappingLine;
		private int missmatches;
		private boolean isSplitLocation;
		private char strandOfSpliceSignal;
		private boolean overlapsKnownJunction;
		private boolean overlapsCleavageSite;
		
		private int start;
		private int end;
		
		public MultiMappingLine(String mappingLine, int missmatches,int start, int end, boolean isSplitLocation, char strandOfSpliceSignal, boolean overlapsKnownJunction, boolean overlapsCleavageSite) {
			this.mappingLine = mappingLine;
			this.missmatches = missmatches;
			this.isSplitLocation = isSplitLocation;
			this.strandOfSpliceSignal = strandOfSpliceSignal;
			this.overlapsKnownJunction = overlapsKnownJunction;
			this.overlapsCleavageSite = overlapsCleavageSite;
			this.start = start;
			this.end = end;
		}
		
		public void update(String mappingLine, int missmatches,int start, int end, boolean isSplitLocation, char strandOfSpliceSignal, boolean overlapsKnownJunction, boolean overlapsCleavageSite) {
			this.mappingLine = mappingLine;
			this.missmatches = missmatches;
			this.isSplitLocation = isSplitLocation;
			this.strandOfSpliceSignal = strandOfSpliceSignal;
			this.overlapsKnownJunction = overlapsKnownJunction;
			this.overlapsCleavageSite = overlapsCleavageSite;
			this.start = start;
			this.end = end;
		}
		
		public String getLine() {
			return this.mappingLine;
		}
		
		public int getMissmatches() {
			return this.missmatches;
		}
		
		public int getStart() {
			return this.start;
		}
		
		public int getEnd() { 
			return this.end;
		}
		
		public boolean isSplitLocation() {
			return this.isSplitLocation;
		}
		
		public char getStrandOfSpliceSignal() {
			return this.strandOfSpliceSignal;
		}
		
		public boolean overlapsKnownJunction() {
			return this.overlapsKnownJunction;
		}
		
		public boolean overlapsCleavageSite() {
			return this.overlapsCleavageSite;
		}
	
		@Override
		public int compareTo(Object o) {
			int missmatchesA = this.missmatches;
			int missmatchesB = ((MultiMappingLine)o).getMissmatches();
			return Integer.valueOf(missmatchesA).compareTo(missmatchesB);
		}
		
	}
	
	
	//expects the best matching alignments,splitted by chr and sorted by start position
	public void resolveOverlappingSpliceSites(String mappingFilePath, int readLength,int maxMissmatches,StringBuilder localReference, int contextOffset, String outputPath) {
		try {
			UnsynchronizedBufferedWriter pw = new UnsynchronizedBufferedWriter(new FileWriter(new File(outputPath)),10240);
			
			
			ArrayList<ExtendedReadLocation> currentSplitLocations = new ArrayList<ExtendedReadLocation>();
			ArrayList<ExtendedReadLocation> locationPool = new ArrayList<ExtendedReadLocation>();
			ExtendedReadLocation tmpLocation;
			
			ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices = new ArrayList<Triplet<Integer,Integer,Integer>>(); 
			ArrayList<ExtendedReadLocation> currentClippedLocationsAtStart = new ArrayList<ExtendedReadLocation>();
			ArrayList<ExtendedReadLocation> currentClippedLocationsAtEnd = new ArrayList<ExtendedReadLocation>();
			ExtendedReadLocation tmpClippedLocation;
			
			
			
			
			String currentLine;
			StringTokenizer st;
			
			String contextId;
			String readId;
			String chr;
			char strand;
			int startA = 0;
			int lastStartClippedAtStart = Integer.MIN_VALUE;
			int lastStartClippedAtEnd = Integer.MIN_VALUE;
			int endA = 0;
			int lastSplitPoint = 0;
			int lastClippingPointAtStart = 0;
			int lastClippingPointAtEnd = 0;
			int startB = 0;
			int endB = 0;
			int upstreamClippingLength;
			int downstreamClippingIndex;
			int downstreamClippingLength;
			int fullAlignmentMismatches;
			int mismatches;
			int overallMappingCount;
			int overallValidPairCount;
			String readSequence;
			char strandOfSpliceSignal;
			boolean overlapsKnownJunction;
			String fullReadClippingLine = null;
			//first element: no polyA tail. second element: reads with a tail
			Pair<ArrayList<ExtendedReadLocation>,ArrayList<ExtendedReadLocation>> splittedOverlappingClippedReads;
			StringBuilder readSequenceBuffer = new StringBuilder();
			
			String polyAtail = "";
			String polyAtailRc = "";
			for(int i = 0; i < this.minPolyALength; i++) {
				polyAtail += "A";
				polyAtailRc += "T";
			}
			
			long prevFilePointer;
			
			
			//if a full alignment has more than 3 mismatches we check if we get a better alignment score with a clipped alignment (only when polyA prediction is enabled)
			if(this.polyA) {
				BufferedReader tmpReader = new BufferedReader(new FileReader(new File(mappingFilePath)));
				String newMappingFilePath = new File(mappingFilePath).getParent() + "/multi_mapping.best.matchings.txt.startposition.sorted.clipped";
				PrintWriter tmpWriter = new PrintWriter(new FileWriter(new File(newMappingFilePath)));
				
				while((currentLine = tmpReader.readLine()) != null) {
					
					//we always print the current line. then we check if we can add for a full read alignment a clipped alignment
					tmpWriter.println(currentLine);
					
				
				
					st = new StringTokenizer(currentLine,"\t");
					//context id
					contextId = st.nextToken();
					//read id
					readId = st.nextToken();
					//chr
					chr = st.nextToken();
					//strand
					strand = st.nextToken().toCharArray()[0];
					
					startA = Integer.valueOf(st.nextToken());
					endA = Integer.valueOf(st.nextToken());
					startB = Integer.valueOf(st.nextToken());
					endB = Integer.valueOf(st.nextToken());
					mismatches = Integer.valueOf(st.nextToken());
					
					if(startB == 0 && endB == 0 && mismatches > 3) {
					
						strandOfSpliceSignal = st.nextToken().charAt(0);
						overlapsKnownJunction = Boolean.valueOf(st.nextToken());
						overallMappingCount = Integer.valueOf(st.nextToken());
						overallValidPairCount = Integer.valueOf(st.nextToken());
						readLength = endA - startA + 1;
						upstreamClippingLength = startB * -1;
						downstreamClippingIndex = endB * -1;
						downstreamClippingLength = readLength - (downstreamClippingIndex  + 1);
					
						
						readSequence = this.string2bitset.decompress(this.read2sequence.get(readId));
						clippingIndices =  getBestLocalAlignmentIndices(readSequence,localReference,contextOffset, readSequenceBuffer, startA, startA + readSequence.length() -1 ,strand, 1, -4, maxMissmatches,this.seedLength,false,false);
						
						for(int j = 0; j < clippingIndices.size();j++) {
							startB = -1 * clippingIndices.get(j).getFirst();
							endB = -1 * clippingIndices.get(j).getSecond();
							upstreamClippingLength = startB * -1;
							downstreamClippingIndex = endB * -1;
							downstreamClippingLength = readLength - (downstreamClippingIndex  + 1);
							fullAlignmentMismatches = mismatches;
							mismatches = clippingIndices.get(j).getThird();
							//fewer mismatches than the full alignment, clipped but not clipped at both ends
							if(mismatches < fullAlignmentMismatches && (upstreamClippingLength >= 3 || (downstreamClippingIndex > 0 && downstreamClippingLength >= 3)) && !(upstreamClippingLength > 0 && downstreamClippingIndex > 0 && downstreamClippingLength > 0)) {
								tmpWriter.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",contextId,readId,chr,strand,startA,endA,startB,endB,mismatches,strandOfSpliceSignal,overlapsKnownJunction,overallMappingCount,overallValidPairCount));
							}
						}
					}
				}
				tmpReader.close();
				tmpWriter.close();
				mappingFilePath = newMappingFilePath;
			}
			

			
			File mappingFile = new File(mappingFilePath);
			Pair<HashMap<Integer,HashMap<Integer,MutableInt>>,HashMap<Integer,HashMap<Integer,MutableInt>>> splitAndClippedCoverages = getSplitCoverages(mappingFile,maxMissmatches);
			HashMap<Integer,HashMap<Integer,MutableInt>> splitCoverages = splitAndClippedCoverages.getFirst();
			HashMap<Integer,Double> splitScores = getSplitScores(splitCoverages,maxMissmatches);
			splitCoverages.clear();
			
			HashMap<Integer,HashMap<Integer,MutableInt>> clippedCoverages = splitAndClippedCoverages.getSecond();
			HashMap<Integer,Double> clippedScores = getSplitScores(clippedCoverages,maxMissmatches);
			clippedCoverages.clear();

			BufferedRandomAccessFile br = new BufferedRandomAccessFile(mappingFile, "r",4096);
			BufferedRandomAccessFile brForRandomSeeks = new BufferedRandomAccessFile(mappingFile, "r",4096);
			prevFilePointer = br.getFilePointer();
			while((currentLine = br.getNextLine()) != null) {
				
				st = new StringTokenizer(currentLine,"\t");
				//context id
				contextId = st.nextToken();
				//read id
				readId = st.nextToken();
				//chr
				chr = st.nextToken();
				//strand
				strand = st.nextToken().toCharArray()[0];
				
				startA = Integer.valueOf(st.nextToken());
				endA = Integer.valueOf(st.nextToken());
				startB = Integer.valueOf(st.nextToken());
				endB = Integer.valueOf(st.nextToken());
				mismatches = Integer.valueOf(st.nextToken());
				strandOfSpliceSignal = st.nextToken().charAt(0);
				overlapsKnownJunction = Boolean.valueOf(st.nextToken());
				overallMappingCount = Integer.valueOf(st.nextToken());
				overallValidPairCount = Integer.valueOf(st.nextToken());
				
				readLength = endA - startA + 1;
				upstreamClippingLength = startB * -1;
				downstreamClippingIndex = endB * -1;
				downstreamClippingLength = readLength - (downstreamClippingIndex  + 1);
				
				//complete read
				if((startB == 0 && endB == 0) ||
				//clipped, but polyA prediction deactivated						
				(!this.polyA && startB <= 0) || 
				//clipped at both ends. not used for a prediction
				(upstreamClippingLength > 0 && downstreamClippingIndex > 0 && downstreamClippingLength > 0)) { 
					pw.write(currentLine);
					//polyA overlap -> true: 1, false: 0
					pw.write("\t0");
					pw.newLine();
					prevFilePointer = br.getFilePointer();
					continue;
				}
				
				//clipping only at the start
				else if(upstreamClippingLength > 0) {
					//if(startA > lastClippingPointAtStart || currentClippedLocationsAtStart.size() > 100000) {
					if(startA > lastStartClippedAtStart + readLength - 1 || currentClippedLocationsAtStart.size() > 100000) {

						if(currentClippedLocationsAtStart.size() < this.minPolyAReadCount) {
							for(ExtendedReadLocation erl : currentClippedLocationsAtStart) {
								printLocation(erl,pw,brForRandomSeeks,false);
							}
						}
						
						else {
							splittedOverlappingClippedReads = seperateNormalAndPolyATailReads(currentClippedLocationsAtStart,this.minPolyAReadCount, this.upperCutoff, this.lowerCutoff,this.maxConsideredClippingLength,true,readSequenceBuffer,polyAtail,polyAtailRc,this.strandedPolyA);
							
							if(!splittedOverlappingClippedReads.getFirst().isEmpty()) {
								for(ExtendedReadLocation erl : splittedOverlappingClippedReads.getFirst()) {
									printLocation(erl,pw,brForRandomSeeks,false);
								}
							}
							
							if(!splittedOverlappingClippedReads.getSecond().isEmpty()) {
								//determine exact cleavage site position
								getAndPrintHighScoringPolyAlocations(splittedOverlappingClippedReads.getSecond(), clippedScores,pw, brForRandomSeeks);
							}
							
						}
						
						
						lastClippingPointAtStart = startA + upstreamClippingLength;
						lastStartClippedAtStart = startA;
						currentClippedLocationsAtStart.clear();
					}
					
					tmpClippedLocation = new ExtendedReadLocation(readId,startA,endA,strand,prevFilePointer);
					
					tmpClippedLocation.setClippingLength(upstreamClippingLength);
					currentClippedLocationsAtStart.add(tmpClippedLocation);
					
				}
				
				//clipping only at the end
				else if(downstreamClippingIndex > 0 && downstreamClippingLength > 0) {
					
					//if(endA - downstreamClippingLength > lastClippingPointAtEnd || currentClippedLocationsAtEnd.size() > 100000) {
					if(startA > lastStartClippedAtEnd + readLength - 1 || currentClippedLocationsAtEnd.size() > 100000) {
						if(currentClippedLocationsAtEnd.size() < this.minPolyAReadCount) {
							for(ExtendedReadLocation erl : currentClippedLocationsAtEnd) {
								printLocation(erl,pw,brForRandomSeeks,false);
							}
						}
						
						else {
							splittedOverlappingClippedReads = seperateNormalAndPolyATailReads(currentClippedLocationsAtEnd,this.minPolyAReadCount, this.upperCutoff, this.lowerCutoff,this.maxConsideredClippingLength, false,readSequenceBuffer,polyAtail,polyAtailRc,this.strandedPolyA);
							if(!splittedOverlappingClippedReads.getFirst().isEmpty()) {
								for(ExtendedReadLocation erl : splittedOverlappingClippedReads.getFirst()) {
									printLocation(erl,pw,brForRandomSeeks,false);
								}
							}
							
							if(!splittedOverlappingClippedReads.getSecond().isEmpty()) {
								//determine exact cleavage site position
								getAndPrintHighScoringPolyAlocations(splittedOverlappingClippedReads.getSecond(), clippedScores,pw, brForRandomSeeks);
							}
						}
						
						
						lastClippingPointAtEnd = endA;
						lastStartClippedAtEnd = startA;
						currentClippedLocationsAtEnd.clear();
					}
					
					tmpClippedLocation = new ExtendedReadLocation(readId,startA,endA,strand,prevFilePointer);
					
					tmpClippedLocation.setClippingLength(downstreamClippingLength);
					currentClippedLocationsAtEnd.add(tmpClippedLocation);
					
				}
				
				//split read
				else {
					readLength = (endA - startA + 1) + (endB - startB + 1);
					
					if(Math.abs(endA - lastSplitPoint) > readLength || currentSplitLocations.size() > 100000) {
						
						if(currentSplitLocations.size() == 1) {
							printLocation(currentSplitLocations.get(0),pw,brForRandomSeeks,false);
						}
						
						else if(currentSplitLocations.size() > 0)
							resolveSpliceSites(currentSplitLocations,splitScores,readLength,pw,brForRandomSeeks);
						
						lastSplitPoint = endA;
						currentSplitLocations.clear();
					}
					
					
					if(currentSplitLocations.size() < locationPool.size()) {
						tmpLocation = locationPool.get(currentSplitLocations.size());
						tmpLocation.updateLocation(readId,endA,startB,strand,prevFilePointer);
					}
					else {
						tmpLocation = new ExtendedReadLocation(readId,endA,startB,strand,prevFilePointer);
						locationPool.add(tmpLocation);
					}
					currentSplitLocations.add(tmpLocation);
				}
			
				
				prevFilePointer = br.getFilePointer();
			}
			
			
			//check the last parsed splits
			if(currentSplitLocations.size() == 1) {
				printLocation(currentSplitLocations.get(0),pw,brForRandomSeeks,false);
			}
			
			else if(currentSplitLocations.size() > 0)
				resolveSpliceSites(currentSplitLocations,splitScores,readLength,pw,brForRandomSeeks);
			
			
			//check the last parsed alignments which are clipped at the start only
			if(currentClippedLocationsAtStart.size() < this.minPolyAReadCount) {
				for(ExtendedReadLocation erl : currentClippedLocationsAtStart) {
					printLocation(erl,pw,brForRandomSeeks,false);
				}
			}
			
			else {
				splittedOverlappingClippedReads = seperateNormalAndPolyATailReads(currentClippedLocationsAtStart,this.minPolyAReadCount, this.upperCutoff, this.lowerCutoff,this.maxConsideredClippingLength,true,readSequenceBuffer,polyAtail,polyAtailRc,this.strandedPolyA);
				if(!splittedOverlappingClippedReads.getFirst().isEmpty()) {
					for(ExtendedReadLocation erl : splittedOverlappingClippedReads.getFirst()) {
						printLocation(erl,pw,brForRandomSeeks,false);
					}
				}
				
				if(!splittedOverlappingClippedReads.getSecond().isEmpty()) {
					//process overlapping clippings here
					getAndPrintHighScoringPolyAlocations(splittedOverlappingClippedReads.getSecond(), clippedScores,pw, brForRandomSeeks);
				}
			}
			
			//check the last parsed alignments which are clipped at the end only					
			if(currentClippedLocationsAtEnd.size() < this.minPolyAReadCount) {
				for(ExtendedReadLocation erl : currentClippedLocationsAtEnd) {
					printLocation(erl,pw,brForRandomSeeks,false);
				}
			}
			
			else {
				splittedOverlappingClippedReads = seperateNormalAndPolyATailReads(currentClippedLocationsAtEnd,this.minPolyAReadCount, this.upperCutoff, this.lowerCutoff,this.maxConsideredClippingLength,false,readSequenceBuffer,polyAtail,polyAtailRc,this.strandedPolyA);
				if(!splittedOverlappingClippedReads.getFirst().isEmpty()) {
					for(ExtendedReadLocation erl : splittedOverlappingClippedReads.getFirst()) {
						printLocation(erl,pw,brForRandomSeeks,false);
					}
				}
				
				if(!splittedOverlappingClippedReads.getSecond().isEmpty()) {
					//process overlapping clippings here
					getAndPrintHighScoringPolyAlocations(splittedOverlappingClippedReads.getSecond(), clippedScores,pw, brForRandomSeeks);
				}
			}
				
			lastClippingPointAtEnd = 0;
			currentClippedLocationsAtEnd.clear();
			lastClippingPointAtStart = 0;
			currentClippedLocationsAtStart.clear();
			currentSplitLocations.clear();
			lastSplitPoint = 0;
			br.close();
			brForRandomSeeks.close();
		
			
			pw.close();
			currentClippedLocationsAtStart.clear();
			currentClippedLocationsAtEnd.clear();
			currentSplitLocations.clear();
			locationPool.clear();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	//first element: reads without polyA-tail
	//second element: reads with polyA-tail
	
	private Pair<ArrayList<ExtendedReadLocation>,ArrayList<ExtendedReadLocation>> seperateNormalAndPolyATailReads(ArrayList<ExtendedReadLocation> readLocations, int minPolyAReadCount, double upperCutoff, double lowerCutoff, int maxConsideredClippingLength, boolean clippedAtStart,StringBuilder readSequenceBuffer, String polyAtail, String polyAtailRc, boolean strandedPolyA) {

		Pair<ArrayList<ExtendedReadLocation>,ArrayList<ExtendedReadLocation>> result = new Pair<ArrayList<ExtendedReadLocation>,ArrayList<ExtendedReadLocation>>();
		ArrayList<ExtendedReadLocation> noTail = new ArrayList<ExtendedReadLocation>();
		ArrayList<ExtendedReadLocation> withTail = new ArrayList<ExtendedReadLocation>();
		result.setFirst(noTail);
		result.setSecond(withTail);
		
		int firstClippingPosition = Integer.MAX_VALUE;
		int lastClippingPosition = Integer.MIN_VALUE;
		int meanClippingLength = 0;
		int smallestAlignmentStart= Integer.MAX_VALUE;
		int largestAlignmentEnd = Integer.MIN_VALUE;
		
		/*
		 * 
		 * check if there are at least minPolyALength consecutive A's or T's in every 
		 * sequence of the overlapping clippings
		 * 
		 */
		String tmpSequence;
		int readLength;
		int aCounter = 0;
		int tCounter = 0;
		boolean foundPolyAWindow = false;
		HashSet<Integer> clippingPositionsWithLongEnoughClippings = new HashSet<Integer>();
		StringBuffer windowRef = new StringBuffer(); 
		
		//actually the ExtendedReadLocation object contains the split coordinates [endA,startB] in the respective fields endA and startB.
		//However, for clipped alignments we save in endA the unclipped read start and in startB the unclipped read end.
		//Thus startB - endA + 1 = readLength holds.
		
		
		Collections.sort(readLocations, Collections.reverseOrder(new ExtendedReadLocationClippingLengthComparator()));
		for(ExtendedReadLocation l : readLocations) {
			
			readLength = l.getStartB() - l.getEndA() + 1;
			tmpSequence = this.string2bitset.decompress(this.read2sequence.get(l.getReadId())).toUpperCase();
			if(l.getStrand() == '-') {
				readSequenceBuffer.setLength(0);
				readSequenceBuffer.append(tmpSequence);
				readSequenceBuffer.reverse();
				for(int i = 0; i < readSequenceBuffer.length(); i++) {
					readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
				}
				tmpSequence = readSequenceBuffer.toString();
			}
			
			if(!clippedAtStart) {
				
				l.setPositionBeforeClippingStarts(l.getEndA() + readLength - l.getClippingLength() - 1);
				tmpSequence = tmpSequence.substring(readLength - l.getClippingLength());
				
				if(tmpSequence.length() > maxConsideredClippingLength)
					tmpSequence = tmpSequence.substring(0,maxConsideredClippingLength);
				
				//check in windows if more than cutoff bases are A's or T's
				foundPolyAWindow = false;
				if(tmpSequence.length() >= this.minPolyALength) {
					for(int i = 0; i <= tmpSequence.length() - polyAtail.length(); i++) {
						windowRef.setLength(0);
						windowRef.append(tmpSequence.substring(i,i + polyAtail.length()));
						aCounter = 0;
						tCounter = 0;
						for(int j = 0; j < windowRef.length(); j++) {
							if(windowRef.charAt(j) == 'A' && (!strandedPolyA || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("1"))))
								aCounter++;
							else if(windowRef.charAt(j) == 'T' && (!strandedPolyA || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("1"))))
								tCounter++;
						}
						
						if((double)aCounter/(double)windowRef.length() >= upperCutoff || (double)tCounter/(double)windowRef.length() >= upperCutoff) {
							foundPolyAWindow = true;
						}
						
						//window with low A or T content. Most possibly a fp hit
						else if((double)aCounter/(double)windowRef.length() <= lowerCutoff && (double)tCounter/(double)windowRef.length() <= lowerCutoff) {
							foundPolyAWindow = false;
							break;
						}
					}
				}
				
				//clipping shorter than the user defined window length. here we require that all bases are A's or T's
				else if(clippingPositionsWithLongEnoughClippings.contains(l.getPositionBeforeClippingStarts()) && l.getClippingLength() > 3) {
					if((polyAtail.contains(tmpSequence) && (!strandedPolyA || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("1")))) || (polyAtailRc.contains(tmpSequence) && (!strandedPolyA || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("1"))))) {
						foundPolyAWindow = true;
					}
				}
				
				
				if(!foundPolyAWindow) {
					noTail.add(l);
				}
				
				else {
					withTail.add(l);
					if(l.getClippingLength() >= this.minPolyALength) {
						
						clippingPositionsWithLongEnoughClippings.add(l.getPositionBeforeClippingStarts());
					}
				}
				
			}
			else {
				l.setPositionBeforeClippingStarts(l.getEndA() + l.getClippingLength());
				tmpSequence = tmpSequence.substring(0,l.getClippingLength());
				
				if(tmpSequence.length() > maxConsideredClippingLength)
					tmpSequence = tmpSequence.substring(tmpSequence.length() - maxConsideredClippingLength);
				
				//check in windows if more than cutoff bases are A's or T's
				foundPolyAWindow = false;
				if(tmpSequence.length() >= this.minPolyALength) {
					for(int i = 0; i <= tmpSequence.length() - polyAtail.length(); i++) {
						windowRef.setLength(0);
						windowRef.append(tmpSequence.substring(i,i + polyAtail.length()));
						aCounter = 0;
						tCounter = 0;
						for(int j = 0; j < windowRef.length(); j++) {
							if(windowRef.charAt(j) == 'A' && (!strandedPolyA || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("1"))))
								aCounter++;
							else if(windowRef.charAt(j) == 'T' && (!strandedPolyA || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("1"))))
								tCounter++;
						}
						
						if((double)aCounter/(double)windowRef.length() >= upperCutoff || (double)tCounter/(double)windowRef.length() >= upperCutoff) {
							foundPolyAWindow = true;
						}
						
						//window with low A or T content. Most possibly a fp hit
						else if((double)aCounter/(double)windowRef.length() <= lowerCutoff && (double)tCounter/(double)windowRef.length() <= lowerCutoff) {
							foundPolyAWindow = false;
							break;
						}
					}
				}
				else if((clippingPositionsWithLongEnoughClippings.contains(l.getPositionBeforeClippingStarts()) && l.getClippingLength() > 3)) {
					if((polyAtail.contains(tmpSequence) && (!strandedPolyA || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("1")))) || (polyAtailRc.contains(tmpSequence) && (!strandedPolyA || (l.getStrand() == '-' && l.getReadId().substring(l.getReadId().length() - 1).equals("2")) || (l.getStrand() == '+' && l.getReadId().substring(l.getReadId().length() - 1).equals("1"))))) {
						foundPolyAWindow = true;
					}
				}
				
				
				
				
				if(!foundPolyAWindow) {
					noTail.add(l);
				}
				
				else {
					withTail.add(l);
					if(l.getClippingLength() >= this.minPolyALength) {
						clippingPositionsWithLongEnoughClippings.add(l.getPositionBeforeClippingStarts());
					}
				}
			}
		}
		
		
		//too few polyA tail reads found -> write out all clippings without polyA-tail flag
		if(withTail.isEmpty() || withTail.size() < minPolyAReadCount) {
			noTail.addAll(withTail);
			withTail.clear();
			return result;
		}
		
		return result;
	}
	
	
	
	@SuppressWarnings("unchecked")
	private void resolveSpliceSites(ArrayList<ExtendedReadLocation> readLocations,HashMap<Integer,Double> splitScores, int readLength, UnsynchronizedBufferedWriter pw, BufferedRandomAccessFile br) throws Exception {
		
		//first check if the second split positions also overlap
		ArrayList<ExtendedReadLocation> unsolvedLocations = new ArrayList<ExtendedReadLocation>();
		Collections.sort(readLocations, new ExtendedReadLocationSecondSplitPositionComparator());
		
		int firstSplitB = readLocations.get(0).getStartB();
		HashMap<Integer,HashSet<Integer>> containedSplits = new HashMap<Integer,HashSet<Integer>>();
		unsolvedLocations.add(readLocations.get(0));
		HashSet<Integer> tmpSet = new HashSet<Integer>();
		tmpSet.add(readLocations.get(0).getStartB());
		containedSplits.put(readLocations.get(0).getEndA(),tmpSet);
		for(int i = 1; i < readLocations.size(); i++) {
			if(readLocations.get(i).getStartB() - firstSplitB < readLength) {
				unsolvedLocations.add(readLocations.get(i));
				if(!containedSplits.containsKey(readLocations.get(i).getEndA()))
					containedSplits.put(readLocations.get(i).getEndA(),new HashSet<Integer>());
				
				containedSplits.get(readLocations.get(i).getEndA()).add(readLocations.get(i).getStartB());
			}
			else {
				
				if(unsolvedLocations.size() == 1) {
					printLocation(unsolvedLocations.get(0),pw,br,false);
				}
				else {
					getAndPrintHighScoringSplitLocations(unsolvedLocations,containedSplits,splitScores,pw,br);
				}
				
				firstSplitB = readLocations.get(i).getStartB();
				containedSplits.clear();
				unsolvedLocations.clear();
				unsolvedLocations.add(readLocations.get(i));
				tmpSet = new HashSet<Integer>();
				tmpSet.add(readLocations.get(i).getStartB());
				containedSplits.put(readLocations.get(i).getEndA(),tmpSet);
			}
		}	
		//process the remaining unsolved locations here
		if(unsolvedLocations.size() == 1)
			printLocation(unsolvedLocations.get(0),pw,br,false);
		
		else if(unsolvedLocations.size() > 1) {
			getAndPrintHighScoringSplitLocations(unsolvedLocations,containedSplits,splitScores,pw,br);
		}
	}
	
	
	private void getAndPrintHighScoringSplitLocations(ArrayList<ExtendedReadLocation> locations,HashMap<Integer,HashSet<Integer>> containedSplits, HashMap<Integer,Double> splitScores, UnsynchronizedBufferedWriter pw, BufferedRandomAccessFile br) throws Exception {
		//first determine best split position
		double currentScore;
		StringBuilder sb = new StringBuilder();
		ArrayList<Split> splits = new ArrayList<Split>();
		SplitScoreComparator scoreComparator = new SplitScoreComparator();
		for(Integer splitPosA : containedSplits.keySet()) {
			for(Integer splitPosB : containedSplits.get(splitPosA)) {
				currentScore = splitScores.get(splitPosA) + splitScores.get(splitPosB);
				Split tmpSplit = new Split(splitPosA,splitPosB);
				tmpSplit.setScore(currentScore);
				splits.add(tmpSplit);
			}
		}
		
		Collections.sort(splits,scoreComparator);
		HashSet<String> printedReadIds = new HashSet<String>();
		for(int i = 0; i < splits.size() && i < 3; i++) {
			for(int j = 0; j < locations.size(); j++) {
				
				if(printedReadIds.contains(locations.get(j).getReadId()))
					continue;
				
				if(splits.get(i).getSplitPositionA() == locations.get(j).getEndA() && splits.get(i).getSplitPositionB() == locations.get(j).getStartB()) {
					br.seek(locations.get(j).getFilePointer());
					pw.write(br.getNextLine() + "\t0");
					pw.newLine();
					printedReadIds.add(locations.get(j).getReadId());
				}
			}
		}
	}
	
	
	private void getAndPrintHighScoringPolyAlocations(ArrayList<ExtendedReadLocation> locations, HashMap<Integer,Double> clippedScores, UnsynchronizedBufferedWriter pw, BufferedRandomAccessFile br) throws Exception {
		//first determine best polyA position
		ArrayList<TranscriptEnd> transcriptEnds = new ArrayList<TranscriptEnd>();
		HashSet<Integer> addedEnds = new HashSet<Integer>();
		for(ExtendedReadLocation l : locations) {
			if(!addedEnds.contains(l.getPositionBeforeClippingStarts())) {
				transcriptEnds.add(new TranscriptEnd(l.getPositionBeforeClippingStarts(),clippedScores.get(l.getPositionBeforeClippingStarts())));
				addedEnds.add(l.getPositionBeforeClippingStarts());
			}
		}
		
		Collections.sort(transcriptEnds);
		HashSet<String> printedReadIds = new HashSet<String>();
		//for(int i = transcriptEnds.size() - 1; i >= 0 && i >= transcriptEnds.size() - 2; i--) {
		for(int i = transcriptEnds.size() - 1; i >= 0; i--) {
			for(int j = 0; j < locations.size(); j++) {
				
				if(printedReadIds.contains(locations.get(j).getReadId()))
					continue;
				
				if(transcriptEnds.get(i).getEnd() == locations.get(j).getPositionBeforeClippingStarts()) {
					br.seek(locations.get(j).getFilePointer());
					pw.write(br.getNextLine() + "\t1");
					pw.newLine();
					printedReadIds.add(locations.get(j).getReadId());
				}
			}
		}
	}
	
	
	private class TranscriptEnd implements Comparable {
		private int end;
		private double score;
		
		public TranscriptEnd(int end, double score) {
			this.end = end;
			this.score = score;
		}

		@Override
		public int compareTo(Object o) {
			// TODO Auto-generated method stub
			double cmp = ((TranscriptEnd) o).getScore();
			return new Double(this.score).compareTo(cmp);
		}
		
		public double getScore() {
			return this.score;
		}
		
		public int getEnd() {
			return this.end;
		}
		
	}
	

	
	private class ExtendedReadLocation {
		private String readId;
		private int endA;
		private int startB;
		
		private char strand;
		
		private int clippingLength;
		private int positionBeforeClippingStarts;
		
		private long filePointer;
		
		public ExtendedReadLocation(String readId, int endA, int startB, char strand, long filePointer) {
			this.readId = readId;
			this.endA = endA;
			this.startB = startB;
			this.strand = strand;
			this.filePointer = filePointer;
		}
		
		public void updateLocation(String readId, int endA, int startB, char strand, long filePointer) {
			this.readId = readId;
			this.endA = endA;
			this.startB = startB;
			this.strand = strand;
			this.filePointer = filePointer;
		}

		public String getReadId() {
			return this.readId;
		}
		
		public int getEndA() {
			return endA;
		}

		public int getStartB() {
			return startB;
		}
		
		
		public char getStrand() {
			return this.strand;
		}
		
		public long getFilePointer() {
			return this.filePointer;
		}
		
		public int getClippingLength() {
			return this.clippingLength;
		}
		
		public void setClippingLength(int clippingLength) {
			this.clippingLength = clippingLength;
		}
		
		public int getPositionBeforeClippingStarts() {
			return this.positionBeforeClippingStarts;
		}
		
		public void setPositionBeforeClippingStarts(int positionBeforeClippingStarts) {
			this.positionBeforeClippingStarts = positionBeforeClippingStarts;
		}
		
	}
	
	private class ExtendedReadLocationSecondSplitPositionComparator implements Comparator<ExtendedReadLocation> {
		
		public int compare(ExtendedReadLocation l1, ExtendedReadLocation l2) {
			int startL1 = l1.getStartB();
			int startL2 = l2.getStartB();
			return ((Integer)startL1).compareTo(startL2);
		}
	}
	
	
	private class ExtendedReadLocationIdComparator implements Comparator<ExtendedReadLocation> {
		
		public int compare(ExtendedReadLocation l1, ExtendedReadLocation l2) {
			String idL1 = l1.getReadId();
			String idL2 = l2.getReadId();
			return (idL1).compareTo(idL2);
		}
	}
	
	private class ExtendedReadLocationClippingLengthComparator implements Comparator<ExtendedReadLocation> {
		
		public int compare(ExtendedReadLocation l1, ExtendedReadLocation l2) {
			int clippingLength1 = l1.getClippingLength();
			int clippingLength2 = l2.getClippingLength();
			return ((Integer)clippingLength1).compareTo(clippingLength2);
		}
	}

	
	
	private class StartPositionComparator implements Comparator {
		public StartPositionComparator() {
			
		}
		@Override
		public int compare(Object o1, Object o2) {
			ReadLocation l1 = (ReadLocation)o1;
			ReadLocation l2 = (ReadLocation)o2;
			return(Integer.valueOf(l1.getStartA()).compareTo(l2.getStartA()));
		}
	}
	
	
	private void printLocation(ExtendedReadLocation location,UnsynchronizedBufferedWriter pw, BufferedRandomAccessFile br, boolean overlapsPolyAtail) throws Exception {
		br.seek(location.getFilePointer());
		pw.write(br.getNextLine());
		
		if(overlapsPolyAtail)
			pw.write("\t1");
		else
			pw.write("\t0");
		
		pw.newLine();
	}
	
	private Pair<HashMap<Integer,HashMap<Integer,MutableInt>>,HashMap<Integer,HashMap<Integer,MutableInt>>> getSplitCoverages(File mappingFile, int maxMissmatches) throws Exception {
		//init the coverage maps
		Pair<HashMap<Integer,HashMap<Integer,MutableInt>>,HashMap<Integer,HashMap<Integer,MutableInt>>> result = new Pair<HashMap<Integer,HashMap<Integer,MutableInt>>,HashMap<Integer,HashMap<Integer,MutableInt>>>();
		HashMap<Integer,HashMap<Integer,MutableInt>> splitCoverages = new HashMap<Integer,HashMap<Integer,MutableInt>>();
		HashMap<Integer,HashMap<Integer,MutableInt>> clippedCoverages = new HashMap<Integer,HashMap<Integer,MutableInt>>();
		result.setFirst(splitCoverages);
		result.setSecond(clippedCoverages);
		
		HashMap<Integer,HashMap<Integer,MutableInt>> completeCoverage = new HashMap<Integer,HashMap<Integer,MutableInt>>();
		
		HashSet<Integer> splitPositions = new HashSet<Integer>();
		HashSet<Integer> clippedPositions = new HashSet<Integer>();
	
		
		BufferedReader br = new BufferedReader(new FileReader(mappingFile));
		String currentLine;
		StringTokenizer st;
		ArrayList<ReadLocation> currentLocations = new ArrayList<ReadLocation>();
		String contextId;
		String readId;
		String chr;
		String strand;
		int startA = 0;
		int endA = 0;
		int lastCoveredPosition = 0;
		int startB = 0;
		int endB = 0;
		int missmatches;
		boolean hasSpliceSignal;
		
		int readLength;
		int upstreamClippingLength;
		int downstreamClippingIndex;
		int downstreamClippingLength;
		
		while(br.ready()) {
			currentLine = br.readLine();
			st = new StringTokenizer(currentLine,"\t");
			contextId = st.nextToken();
			readId = st.nextToken();
			chr = st.nextToken();
			strand = st.nextToken();
			startA = Integer.valueOf(st.nextToken());
			endA = Integer.valueOf(st.nextToken());
			startB = Integer.valueOf(st.nextToken());
			endB = Integer.valueOf(st.nextToken());
			readLength = endA - startA + 1;
			missmatches = Integer.valueOf(st.nextToken());
			hasSpliceSignal = Boolean.valueOf(st.nextToken());
			
			
			//adding the start and end of full/split reads of the first read part
			if(startB >= 0 && endB >= 0) {
				if(completeCoverage.containsKey(startA)) {
					if(completeCoverage.get(startA).containsKey(missmatches)) {
						completeCoverage.get(startA).get(missmatches).increment();
					}
					else {
						completeCoverage.get(startA).put(missmatches, new MutableInt(1));
					}
				}
				else {
					HashMap<Integer,MutableInt> tmpMap = new HashMap<Integer,MutableInt>();
					tmpMap.put(missmatches, new MutableInt(1));
					completeCoverage.put(startA, tmpMap);
				}
				
				
				//adding the end
				if(completeCoverage.containsKey(endA)) {
					if(completeCoverage.get(endA).containsKey(missmatches)) {
						completeCoverage.get(endA).get(missmatches).increment();
					}
					else {
						completeCoverage.get(endA).put(missmatches, new MutableInt(1));
					}
				}
				else {
					HashMap<Integer,MutableInt> tmpMap = new HashMap<Integer,MutableInt>();
					tmpMap.put(missmatches, new MutableInt(1));
					completeCoverage.put(endA, tmpMap);
				}
			}
			
			//found a split read -> also adding start of second read part
			if(startB > 0) {
				if(completeCoverage.containsKey(startB)) {
					if(completeCoverage.get(startB).containsKey(missmatches)) {
						completeCoverage.get(startB).get(missmatches).increment();
					}
					else {
						completeCoverage.get(startB).put(missmatches, new MutableInt(1));
					}
				}
				else {
					HashMap<Integer,MutableInt> tmpMap = new HashMap<Integer,MutableInt>();
					tmpMap.put(missmatches, new MutableInt(1));
					completeCoverage.put(startB, tmpMap);
				}
				
				
				//add split positions to hashset
				splitPositions.add(endA);
				splitPositions.add(startB);
				
			}
			
			//clipping only at the start (startB and endB give the indices where the clipped alignment starts)
			else if(startB < 0 && (-1) * endB == readLength - 1) {
				upstreamClippingLength = startB * -1;
				
				if(completeCoverage.containsKey(startA + upstreamClippingLength)) {
					if(completeCoverage.get(startA + upstreamClippingLength).containsKey(missmatches)) {
						completeCoverage.get(startA + upstreamClippingLength).get(missmatches).increment();
					}
					else {
						completeCoverage.get(startA + upstreamClippingLength).put(missmatches, new MutableInt(1));
					}
				}
				else {
					HashMap<Integer,MutableInt> tmpMap = new HashMap<Integer,MutableInt>();
					tmpMap.put(missmatches, new MutableInt(1));
					completeCoverage.put(startA + upstreamClippingLength, tmpMap);
				}
				
				
				clippedPositions.add(startA + upstreamClippingLength);
			}
			
			
			//clipping only at the end
			else if(startB == 0 && endB < 0) {
				downstreamClippingIndex = endB * -1;
				downstreamClippingLength = readLength - (downstreamClippingIndex  + 1);
				
				if(completeCoverage.containsKey(startA + readLength - downstreamClippingLength - 1)) {
					if(completeCoverage.get(startA + readLength - downstreamClippingLength - 1).containsKey(missmatches)) {
						completeCoverage.get(startA + readLength - downstreamClippingLength - 1).get(missmatches).increment();
					}
					else {
						completeCoverage.get(startA + readLength - downstreamClippingLength - 1).put(missmatches, new MutableInt(1));
					}
				}
				else {
					HashMap<Integer,MutableInt> tmpMap = new HashMap<Integer,MutableInt>();
					tmpMap.put(missmatches, new MutableInt(1));
					completeCoverage.put(startA + readLength - downstreamClippingLength - 1, tmpMap);
				}
				
				clippedPositions.add(startA + readLength - downstreamClippingLength - 1);
			}
			
		}
		br.close();
		
		
		//now fill the split coverages
		for(int i : splitPositions) {
			splitCoverages.put(i, completeCoverage.get(i));
		}
		
		//and the clipped coverages
		for(int i : clippedPositions) {
			clippedCoverages.put(i, completeCoverage.get(i));
		}
		
		
		return result;
	}
	
	private HashMap<Integer,Double> getSplitScores(HashMap<Integer,HashMap<Integer,MutableInt>> splitCoverages, int maxMissmatches) {
		HashMap<Integer,Double> splitScores = new HashMap<Integer,Double>();
		HashMap<Integer,MutableInt> currentCoverageMap;
		double currentScore;
		
		for(int splitPosition : splitCoverages.keySet()) {
			currentCoverageMap = splitCoverages.get(splitPosition);
			currentScore = 0;
			for(int i = 0; i <= maxMissmatches; i++) {
				if(currentCoverageMap.containsKey(i)) {
					currentScore += (Math.pow(0.3, i) * currentCoverageMap.get(i).intValue());
				}
			}
			splitScores.put(splitPosition,currentScore);
			
		}
		return splitScores;
	}
	
	
	public void buildAllMultiSplitCombinations(String multiMappingFilePath,String multiMappingOutputPath, boolean pairedend,StringBuilder localReference,int contextOffset) {
		HashMap<Integer,IntervalTree<ReadLocation>> currentInterval2location = null;
		try {
			
			/**
			 * hash already detected multi splits
			 * 
			 */
			HashMap<String,HashSet<String>> id2splitKeys = new HashMap<String,HashSet<String>>();
			BufferedReader br;
			String currentLine;
			String[] splittedLine;
			String[] splittedCoordinates;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern commaPattern = Pattern.compile(",");
			if(new File(multiMappingOutputPath).exists()) {
				br = new BufferedReader(new FileReader(new File(multiMappingOutputPath)));
				
				
				while((currentLine = br.readLine()) != null) {
					
					splittedLine = tabPattern.split(currentLine);
					splittedCoordinates = commaPattern.split(splittedLine[4]);
					
					if(splittedCoordinates.length <= 4)
						continue;
					
					if(!id2splitKeys.containsKey(splittedLine[1]))
						id2splitKeys.put(splittedLine[1],new HashSet<String>());
					
					id2splitKeys.get(splittedLine[1]).add(splittedLine[4]);
				}
				br.close();
			}
			
			
			br = new BufferedReader(new FileReader(new File(multiMappingFilePath)));
			BufferedWriter multiMappingWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(multiMappingOutputPath),true)), 1024 * 1024 * 10);
			
			ArrayList<Integer> mappingCountsA = new ArrayList<Integer>();
			ArrayList<Integer> mappingCountsB = new ArrayList<Integer>();
			int mappingCount;
			int validPairsCount;
			
			Pair<HashMap<Integer,IntervalTree<ReadLocation>>,HashMap<Integer,IntervalTree<ReadLocation>>> interval2locationPair = new Pair<HashMap<Integer,IntervalTree<ReadLocation>>,HashMap<Integer,IntervalTree<ReadLocation>>>();
			
			interval2locationPair.setFirst(new HashMap<Integer,IntervalTree<ReadLocation>>());
			interval2locationPair.setSecond(new HashMap<Integer,IntervalTree<ReadLocation>>());
			IntervalTree<ReadLocation> currenTree;
			
			
			Pattern doublePointPattern = Pattern.compile("::");
			String readId;
			String prevReadId = "";
			String realReadId = "";
			String tmpReadId;
			String readBaseId = null;
			String prevReadBaseId = "";
			String chr;
			String prevChr = "";
			String[] splittedId;
			int intervalIndex;
			int maxIntervalIndex = Integer.MIN_VALUE;
			ReadLocation tmpLocation;
			ReadLocation startLocation;
			String realReadSequence = "";
			String realReadSequenceRevComp = "";
			String segmentReadSequence;
			ArrayList<Pair<Integer,Integer>> tmpList;
			ArrayList<ReadLocation> currentCombinations = new ArrayList<ReadLocation>();
			ArrayList<ReadLocation> currentClippings = new ArrayList<ReadLocation>();
			ArrayList<String> segmentReadSequences = new ArrayList<String>();
			int startA;
			int endA;
			int startB;
			int endB;
			int upstreamClippingIndex;
			int downstreamClippingIndex;
			int insertionLength;
			int matchLength;
			
			StringBuilder readSequenceBuffer = new StringBuilder();
			
			while((currentLine = br.readLine()) != null) {
				
				splittedLine = tabPattern.split(currentLine);
				readId = splittedLine[1];
				if(readId.contains("MSC::")) {
					splittedId = doublePointPattern.split(readId);
					intervalIndex = Integer.valueOf(splittedId[3]);
					chr = splittedId[2];
					readBaseId = splittedId[0];
					if(pairedend)
						readBaseId = splittedId[0].substring(0,splittedId[0].length() - 2);
					
					
					startA = Integer.valueOf(splittedLine[4]);
					endA = Integer.valueOf(splittedLine[5]);
					startB = Integer.valueOf(splittedLine[6]);
					endB = Integer.valueOf(splittedLine[7]);
					
					tmpList = new ArrayList<Pair<Integer,Integer>>();
					tmpList.add(new Pair<Integer,Integer>(startA,endA));
					tmpList.add(new Pair<Integer,Integer>(startB,endB));
					
					tmpLocation = new ReadLocation(splittedLine[2], splittedLine[3].charAt(0),'S',tmpList, Integer.valueOf(splittedLine[8]));
					tmpLocation.setReadId(readId);
					tmpLocation.setStrandOfSpliceSignal(splittedLine[9].charAt(0));
					tmpLocation.setOverlapsKnownJunction(Boolean.valueOf(splittedLine[10]));
					mappingCount = Integer.valueOf(splittedLine[11]);
					validPairsCount = Integer.valueOf(splittedLine[12]);
					
					
					if(readId.equals(prevReadId) || (readBaseId.equals(prevReadBaseId) && chr.equals(prevChr))) {
						if(!pairedend || readId.substring(readId.length() - 2).equals("/1")) {
							currentInterval2location = interval2locationPair.getFirst();
							
							if(mappingCountsA.isEmpty()) {
								mappingCountsA.add(mappingCount);
								mappingCountsA.add(validPairsCount);
							}
						}
						
						else {
							currentInterval2location = interval2locationPair.getSecond();
							
							if(mappingCountsB.isEmpty()) {
								mappingCountsB.add(mappingCount);
								mappingCountsB.add(validPairsCount);
							}
						}
						
						
						if(!currentInterval2location.containsKey(intervalIndex))
							currentInterval2location.put(intervalIndex, new IntervalTree<ReadLocation>());
						
						try {
							currentInterval2location.get(intervalIndex).add(tmpLocation);
						}
						catch(Exception e) {
							e.printStackTrace();
						}
						maxIntervalIndex = Math.max(intervalIndex, maxIntervalIndex);
					}
					
					else {
						
						if(!interval2locationPair.getFirst().isEmpty()) {
							tmpReadId = interval2locationPair.getFirst().values().iterator().next().first().getReadId();
							realReadId = doublePointPattern.split(tmpReadId)[0];
							realReadSequence = this.string2bitset.decompress(this.read2sequence.get(realReadId));
							realReadSequenceRevComp = getRevComplement(realReadSequence, readSequenceBuffer);
							
							
							
							currentCombinations.clear();
							
							if(interval2locationPair.getFirst().size() < 30) {
								for(int i = 0; i < maxIntervalIndex; i+=2) {
									if(interval2locationPair.getFirst().containsKey(i)) {
										Iterator<ReadLocation> locationIt = interval2locationPair.getFirst().get(i).iterator();
										while(locationIt.hasNext()) {
											startLocation = locationIt.next();
											
											if(startLocation.getCoordinates().get(0).getFirst() < 0) {
													continue;
											}
											
											segmentReadSequences.clear();
											segmentReadSequences.add(this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId())));
											determineCombinations(startLocation, i+1, maxIntervalIndex, interval2locationPair.getFirst(),currentCombinations,(startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1) + (startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1),realReadSequence,realReadSequenceRevComp,readSequenceBuffer,segmentReadSequences,localReference,contextOffset);
										}
									}
								}
								
								//in case we did not find any possible combination try to find an alignment overlapping a single intron
								if(currentCombinations.isEmpty()) {
									currentClippings.clear();
									for(int i = 0; i <= maxIntervalIndex; i++) {
										if(interval2locationPair.getFirst().containsKey(i)) {
											Iterator<ReadLocation> locationIt = interval2locationPair.getFirst().get(i).iterator();
											while(locationIt.hasNext()) {
												startLocation = locationIt.next();
												segmentReadSequences.clear();
												segmentReadSequence = this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId()));
												segmentReadSequences.add(segmentReadSequence);
												generateNonOverlappingCoordintes(startLocation,segmentReadSequences,realReadSequence);
												if(startLocation.getCoordinates() != null) {
													int mismatches = getMismatchCount(startLocation,realReadSequence,localReference,contextOffset);
													if(mismatches <= this.maxMissmatches) {
														startLocation.setMismatches(mismatches);
														currentCombinations.add(startLocation);
													}
													else if(this.clipping && mismatches != Integer.MAX_VALUE) {
														int firstSegmentLength = startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1; 
														ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices =  getBestLocalAlignmentIndices((startLocation.getStrand() == '+')?realReadSequence.substring(0,firstSegmentLength):realReadSequenceRevComp.substring(0,firstSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(0).getFirst(), startLocation.getCoordinates().get(0).getFirst() + firstSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,true,false);
														if(!clippingIndices.isEmpty()) {
															upstreamClippingIndex = clippingIndices.get(0).getFirst();
															int secondSegmentLength = startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1;
															insertionLength = 0;
															//insertion found
															if(startLocation.getCoordinates().get(1).getFirst() <= startLocation.getCoordinates().get(0).getSecond()) {
																insertionLength = (-1 *(startLocation.getCoordinates().get(1).getFirst() - startLocation.getCoordinates().get(0).getSecond() - 1));
															}
															mismatches = clippingIndices.get(0).getThird();
															clippingIndices = getBestLocalAlignmentIndices((startLocation.getStrand() == '+')? realReadSequence.substring(realReadSequence.length() - secondSegmentLength):realReadSequenceRevComp.substring(realReadSequence.length() - secondSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(1).getFirst(), startLocation.getCoordinates().get(1).getFirst() + secondSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,false,true);
															if(!clippingIndices.isEmpty()) {
																downstreamClippingIndex = firstSegmentLength + insertionLength + clippingIndices.get(0).getSecond();
																downstreamClippingIndex = Math.min(realReadSequence.length() - 1, downstreamClippingIndex);
																mismatches += clippingIndices.get(0).getThird();
																if(mismatches <= this.maxMissmatches) {
																	startLocation.addCoordinate(new Pair<Integer,Integer>( -1 * upstreamClippingIndex, -1 * downstreamClippingIndex));
																	startLocation.setMismatches(mismatches);
																	currentClippings.add(startLocation);
																}
															}
														}														
													}
												}
											}
										}
									}
								}
								
								
								//still no combination found -> add best clipping possibilities
								if(this.clipping && currentCombinations.isEmpty() && !currentClippings.isEmpty()) {
									currentCombinations.addAll(currentClippings);
								}
								
								writeBestMappingLocations(currentCombinations,id2splitKeys,splittedLine[0],realReadId,mappingCountsA,multiMappingWriter);
							}
							
						}
						if(!interval2locationPair.getSecond().isEmpty()) {
							tmpReadId = interval2locationPair.getSecond().values().iterator().next().first().getReadId();
							realReadId = doublePointPattern.split(tmpReadId)[0];
							realReadSequence = this.string2bitset.decompress(this.read2sequence.get(realReadId));
							realReadSequenceRevComp = getRevComplement(realReadSequence, readSequenceBuffer);
							currentCombinations.clear();
							
							if(interval2locationPair.getSecond().size() < 30) {
								for(int i = 0; i < maxIntervalIndex; i+=2) {
									if(interval2locationPair.getSecond().containsKey(i)) {
										Iterator<ReadLocation> locationIt = interval2locationPair.getSecond().get(i).iterator();
										while(locationIt.hasNext()) {
											startLocation = locationIt.next();
											if(startLocation.getCoordinates().get(0).getFirst() < 0) {
												continue;
											}
											segmentReadSequences.clear();
											segmentReadSequences.add(this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId())));
											determineCombinations(startLocation, i+1, maxIntervalIndex, interval2locationPair.getSecond(),currentCombinations,(startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1) + (startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1),realReadSequence,realReadSequenceRevComp,readSequenceBuffer,segmentReadSequences,localReference,contextOffset);
										}
									}
								}
								
								//in case we did not find any possible combination try to find an alignment overlapping a single intron
								if(currentCombinations.isEmpty()) {
									currentClippings.clear();
									for(int i = 0; i <= maxIntervalIndex; i++) {
										if(interval2locationPair.getSecond().containsKey(i)) {
											Iterator<ReadLocation> locationIt = interval2locationPair.getSecond().get(i).iterator();
											while(locationIt.hasNext()) {
												startLocation = locationIt.next();
												segmentReadSequences.clear();
												segmentReadSequence = this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId()));
												segmentReadSequences.add(segmentReadSequence);
												generateNonOverlappingCoordintes(startLocation,segmentReadSequences,realReadSequence);
												if(startLocation.getCoordinates() != null) {
													int mismatches = getMismatchCount(startLocation,realReadSequence,localReference,contextOffset);
													if(mismatches <= this.maxMissmatches) {
														startLocation.setMismatches(mismatches);
														currentCombinations.add(startLocation);
													}
													else if(this.clipping && mismatches != Integer.MAX_VALUE) {
														int firstSegmentLength = startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1; 
														ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices =  getBestLocalAlignmentIndices((startLocation.getStrand() == '+')?realReadSequence.substring(0,firstSegmentLength):realReadSequenceRevComp.substring(0,firstSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(0).getFirst(), startLocation.getCoordinates().get(0).getFirst() + firstSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,true,false);
														if(!clippingIndices.isEmpty()) {
															upstreamClippingIndex = clippingIndices.get(0).getFirst();
															int secondSegmentLength = startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1;
															insertionLength = 0;
															//insertion found
															if(startLocation.getCoordinates().get(1).getFirst() <= startLocation.getCoordinates().get(0).getSecond()) {
																insertionLength = (-1 *(startLocation.getCoordinates().get(1).getFirst() - startLocation.getCoordinates().get(0).getSecond() - 1));
															}
															mismatches = clippingIndices.get(0).getThird();
															clippingIndices = getBestLocalAlignmentIndices((startLocation.getStrand() == '+')? realReadSequence.substring(realReadSequence.length() - secondSegmentLength):realReadSequenceRevComp.substring(realReadSequence.length() - secondSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(1).getFirst(), startLocation.getCoordinates().get(1).getFirst() + secondSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,false,true);
															if(!clippingIndices.isEmpty()) {
																downstreamClippingIndex = firstSegmentLength + insertionLength + clippingIndices.get(0).getSecond();
																downstreamClippingIndex = Math.min(realReadSequence.length() - 1, downstreamClippingIndex);
																mismatches += clippingIndices.get(0).getThird();
																if(mismatches <= this.maxMissmatches) {
																	startLocation.addCoordinate(new Pair<Integer,Integer>( -1 * upstreamClippingIndex, -1 * downstreamClippingIndex));
																	startLocation.setMismatches(mismatches);
																	currentClippings.add(startLocation);
																}
															}
														}													
													}
												}
											}
										}
									}
								}
								
								//still no combination found -> add best clipping possibilities
								if(this.clipping && currentCombinations.isEmpty() && !currentClippings.isEmpty()) {
									currentCombinations.addAll(currentClippings);
								}
								
								writeBestMappingLocations(currentCombinations,id2splitKeys,splittedLine[0],realReadId,mappingCountsB,multiMappingWriter);
							}
							
						}
							
						interval2locationPair.getFirst().clear();
						interval2locationPair.getSecond().clear();
						mappingCountsA.clear();
						mappingCountsB.clear();
						prevReadId = readId;
						prevReadBaseId = readBaseId;
						prevChr = chr;
						maxIntervalIndex = intervalIndex;
						
						if(!pairedend || readId.substring(readId.length() - 2).equals("/1")) {
							currentInterval2location = interval2locationPair.getFirst();
							mappingCountsA.add(mappingCount);
							mappingCountsA.add(validPairsCount);
						}
						
						else {
							currentInterval2location = interval2locationPair.getSecond();
							mappingCountsB.add(mappingCount);
							mappingCountsB.add(validPairsCount);
						}
						
						
						if(!currentInterval2location.containsKey(intervalIndex))
							currentInterval2location.put(intervalIndex, new IntervalTree<ReadLocation>());
						
						try {
							currentInterval2location.get(intervalIndex).add(tmpLocation);
						}
						catch(Exception e) {
							e.printStackTrace();
						}
						
					}
				}
				
				
				else {
					
					if(!interval2locationPair.getFirst().isEmpty()) {
						tmpReadId = interval2locationPair.getFirst().values().iterator().next().first().getReadId();
						realReadId = doublePointPattern.split(tmpReadId)[0];
						realReadSequence = this.string2bitset.decompress(this.read2sequence.get(realReadId));
						realReadSequenceRevComp = getRevComplement(realReadSequence, readSequenceBuffer);
						currentCombinations.clear();
						
						if(interval2locationPair.getFirst().size() < 30) {
							for(int i = 0; i < maxIntervalIndex; i+=2) {
								if(interval2locationPair.getFirst().containsKey(i)) {
									Iterator<ReadLocation> locationIt = interval2locationPair.getFirst().get(i).iterator();
									while(locationIt.hasNext()) {
										startLocation = locationIt.next();
										if(startLocation.getCoordinates().get(0).getFirst() < 0) {
											continue;
										}
										
										segmentReadSequences.clear();
										segmentReadSequences.add(this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId())));
										determineCombinations(startLocation, i+1, maxIntervalIndex, interval2locationPair.getFirst(),currentCombinations,(startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1) + (startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1),realReadSequence,realReadSequenceRevComp,readSequenceBuffer,segmentReadSequences,localReference,contextOffset);
									}
								}
							}
							
							//in case we did not find any possible combination try to find an alignment overlapping a single intron
							if(currentCombinations.isEmpty()) {
								currentClippings.clear();
								for(int i = 0; i <= maxIntervalIndex; i++) {
									if(interval2locationPair.getFirst().containsKey(i)) {
										Iterator<ReadLocation> locationIt = interval2locationPair.getFirst().get(i).iterator();
										while(locationIt.hasNext()) {
											startLocation = locationIt.next();
											segmentReadSequences.clear();
											segmentReadSequence = this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId()));
											segmentReadSequences.add(segmentReadSequence);
											generateNonOverlappingCoordintes(startLocation,segmentReadSequences,realReadSequence);
											if(startLocation.getCoordinates() != null) {
												int mismatches = getMismatchCount(startLocation,realReadSequence,localReference,contextOffset);
												if(mismatches <= this.maxMissmatches) {
													startLocation.setMismatches(mismatches);
													currentCombinations.add(startLocation);
												}
												else if(this.clipping && mismatches != Integer.MAX_VALUE) {
													int firstSegmentLength = startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1; 
													ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices =  getBestLocalAlignmentIndices((startLocation.getStrand() == '+')?realReadSequence.substring(0,firstSegmentLength):realReadSequenceRevComp.substring(0,firstSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(0).getFirst(), startLocation.getCoordinates().get(0).getFirst() + firstSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,true,false);
													if(!clippingIndices.isEmpty()) {
														upstreamClippingIndex = clippingIndices.get(0).getFirst();
														int secondSegmentLength = startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1;
														insertionLength = 0;
														//insertion found
														if(startLocation.getCoordinates().get(1).getFirst() <= startLocation.getCoordinates().get(0).getSecond()) {
															insertionLength = (-1 *(startLocation.getCoordinates().get(1).getFirst() - startLocation.getCoordinates().get(0).getSecond() - 1));
														}
														mismatches = clippingIndices.get(0).getThird();
														clippingIndices = getBestLocalAlignmentIndices((startLocation.getStrand() == '+')? realReadSequence.substring(realReadSequence.length() - secondSegmentLength):realReadSequenceRevComp.substring(realReadSequence.length() - secondSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(1).getFirst(), startLocation.getCoordinates().get(1).getFirst() + secondSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,false,true);
														if(!clippingIndices.isEmpty()) {
															downstreamClippingIndex = firstSegmentLength + insertionLength + clippingIndices.get(0).getSecond();
															downstreamClippingIndex = Math.min(realReadSequence.length() - 1, downstreamClippingIndex);
															mismatches += clippingIndices.get(0).getThird();
															if(mismatches <= this.maxMissmatches) {
																startLocation.addCoordinate(new Pair<Integer,Integer>( -1 * upstreamClippingIndex, -1 * downstreamClippingIndex));
																startLocation.setMismatches(mismatches);
																currentClippings.add(startLocation);
															}
														}
													}													
												}
											}
										}
									}
								}
							}
							
							//still no combination found -> add best clipping possibilities
							if(this.clipping && currentCombinations.isEmpty() && !currentClippings.isEmpty()) {
								currentCombinations.addAll(currentClippings);
							}
							
							writeBestMappingLocations(currentCombinations,id2splitKeys,splittedLine[0],realReadId,mappingCountsA,multiMappingWriter);
						
						}
					}
					
					if(!interval2locationPair.getSecond().isEmpty()) {
						tmpReadId = interval2locationPair.getSecond().values().iterator().next().first().getReadId();
						realReadId = doublePointPattern.split(tmpReadId)[0];
						realReadSequence = this.string2bitset.decompress(this.read2sequence.get(realReadId));
						realReadSequenceRevComp = getRevComplement(realReadSequence, readSequenceBuffer);
						currentCombinations.clear();
						
						if(interval2locationPair.getSecond().size() < 30) {
							for(int i = 0; i < maxIntervalIndex; i+=2) {
								if(interval2locationPair.getSecond().containsKey(i)) {
									Iterator<ReadLocation> locationIt = interval2locationPair.getSecond().get(i).iterator();
									while(locationIt.hasNext()) {
										startLocation = locationIt.next();
										if(startLocation.getCoordinates().get(0).getFirst() < 0) {
											continue;
										}
										segmentReadSequences.clear();
										segmentReadSequences.add(this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId())));
										determineCombinations(startLocation, i+1, maxIntervalIndex, interval2locationPair.getSecond(),currentCombinations,(startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1) + (startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1),realReadSequence,realReadSequenceRevComp,readSequenceBuffer,segmentReadSequences,localReference,contextOffset);
									}
								}
							}
							
							//in case we did not find any possible combination try to find an alignment overlapping a single intron
							if(currentCombinations.isEmpty()) {
								currentClippings.clear();
								for(int i = 0; i <= maxIntervalIndex; i++) {
									if(interval2locationPair.getSecond().containsKey(i)) {
										Iterator<ReadLocation> locationIt = interval2locationPair.getSecond().get(i).iterator();
										while(locationIt.hasNext()) {
											startLocation = locationIt.next();
											segmentReadSequences.clear();
											segmentReadSequence = this.string2bitset.decompress(this.read2sequence.get(startLocation.getReadId()));
											segmentReadSequences.add(segmentReadSequence);
											generateNonOverlappingCoordintes(startLocation,segmentReadSequences,realReadSequence);
											if(startLocation.getCoordinates() != null) {
												int mismatches = getMismatchCount(startLocation,realReadSequence,localReference,contextOffset);
												if(mismatches <= this.maxMissmatches) {
													startLocation.setMismatches(mismatches);
													currentCombinations.add(startLocation);
												}
												else if(this.clipping && mismatches != Integer.MAX_VALUE) {
													int firstSegmentLength = startLocation.getCoordinates().get(0).getSecond() - startLocation.getCoordinates().get(0).getFirst() + 1; 
													ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices =  getBestLocalAlignmentIndices((startLocation.getStrand() == '+')?realReadSequence.substring(0,firstSegmentLength):realReadSequenceRevComp.substring(0,firstSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(0).getFirst(), startLocation.getCoordinates().get(0).getFirst() + firstSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,true,false);
													if(!clippingIndices.isEmpty()) {
														upstreamClippingIndex = clippingIndices.get(0).getFirst();
														int secondSegmentLength = startLocation.getCoordinates().get(1).getSecond() - startLocation.getCoordinates().get(1).getFirst() + 1;
														insertionLength = 0;
														//insertion found
														if(startLocation.getCoordinates().get(1).getFirst() <= startLocation.getCoordinates().get(0).getSecond()) {
															insertionLength = (-1 *(startLocation.getCoordinates().get(1).getFirst() - startLocation.getCoordinates().get(0).getSecond() - 1));
														}
														mismatches = clippingIndices.get(0).getThird();
														clippingIndices = getBestLocalAlignmentIndices((startLocation.getStrand() == '+')? realReadSequence.substring(realReadSequence.length() - secondSegmentLength):realReadSequenceRevComp.substring(realReadSequence.length() - secondSegmentLength),localReference,contextOffset, readSequenceBuffer, startLocation.getCoordinates().get(1).getFirst(), startLocation.getCoordinates().get(1).getFirst() + secondSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,false,true);
														if(!clippingIndices.isEmpty()) {
															downstreamClippingIndex = firstSegmentLength + insertionLength + clippingIndices.get(0).getSecond();
															downstreamClippingIndex = Math.min(realReadSequence.length() - 1, downstreamClippingIndex);
															mismatches += clippingIndices.get(0).getThird();
															
															if(mismatches <= this.maxMissmatches) {
																startLocation.addCoordinate(new Pair<Integer,Integer>( -1 * upstreamClippingIndex, -1 * downstreamClippingIndex));
																startLocation.setMismatches(mismatches);
																currentClippings.add(startLocation);
															}
														}
													}														
												}
											}
										}
									}
								}
							}
							
							//still no combination found -> add best clipping possibilities
							if(this.clipping && currentCombinations.isEmpty() && !currentClippings.isEmpty()) {
								currentCombinations.addAll(currentClippings);
							}
							
							writeBestMappingLocations(currentCombinations,id2splitKeys,splittedLine[0],realReadId,mappingCountsB,multiMappingWriter);
						}
					}
					
					/**
					 * format switch
					 * 
					 * old: context_id	read_id	chr	strand	startA	endA	startB	endB	mismatches	known_splicesignal	known_junction	mappingCount	validPairsCount	overlapsPolyATail
					 * new:	context_id	read_id	chr	strand	startA,endA,startB,endB,...,startN,endN	mismatches	known_splicesignal	known_junction	mappingCount	validPairsCount	overlapsPolyATail
					 */
					
					//multiMappingWriter.write(String.format("%s\t%s\t%s\t%s\t%s,%s,%s,%s\t%s\t%s\t%s\t%s\t%s\t%s\n", splittedLine[0],splittedLine[1],splittedLine[2],splittedLine[3],splittedLine[4],splittedLine[5],splittedLine[6],splittedLine[7],splittedLine[8],splittedLine[9],splittedLine[10],splittedLine[11],splittedLine[12],splittedLine[13]));
					multiMappingWriter.write(splittedLine[0].concat("\t").concat(splittedLine[1]).concat("\t").concat(splittedLine[2]).concat("\t").concat(splittedLine[3]).concat("\t").concat(splittedLine[4]).concat(",").concat(splittedLine[5]).concat(",").concat(splittedLine[6]).concat(",").concat(splittedLine[7]).concat("\t").concat(splittedLine[8]).concat("\t").concat(splittedLine[9]).concat("\t").concat(splittedLine[10]).concat("\t").concat(splittedLine[11]).concat("\t").concat(splittedLine[12]).concat("\t").concat(splittedLine[13]).concat("\n"));
					interval2locationPair.getFirst().clear();
					interval2locationPair.getSecond().clear();
					
					prevReadId = "";
					prevReadBaseId = "";
					prevChr = "";
					maxIntervalIndex = Integer.MIN_VALUE;
				}
			}
			
			multiMappingWriter.close();
			br.close();
		}
		
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	private String getRevComplement(String readSequence, StringBuilder readSequenceBuffer) {
		readSequenceBuffer.setLength(0);
		readSequenceBuffer.append(readSequence);
		readSequenceBuffer.reverse();
		for(int i = 0; i < readSequenceBuffer.length(); i++) {
			readSequenceBuffer.setCharAt(i, substitute(readSequenceBuffer.charAt(i)));
		}
		return readSequenceBuffer.toString();
	}
	
	
	private void writeBestMappingLocations(ArrayList<ReadLocation> currentCombinations, HashMap<String,HashSet<String>> id2splitKeys, String contextId, String tmpReadId,ArrayList<Integer> mappingCounts, BufferedWriter multiMappingWriter) throws Exception {
		//determine min mismatch value
		int minMismatchValue = Integer.MAX_VALUE;
		for(ReadLocation location : currentCombinations) {
			if(location.getMismatches() < minMismatchValue)
				minMismatchValue = location.getMismatches();
		}
		
		StringBuffer sb = new StringBuffer();
		HashSet<String> writtenSplitCandidates = new HashSet<String>();
		for(ReadLocation location : currentCombinations) {
			if(location.getMismatches() - minMismatchValue <= this.maxMissmatchDifference) {
				
				sb.setLength(0);
				sb.append(location.getCoordinates().get(0).getFirst()).append(",").append(location.getCoordinates().get(0).getSecond());
				for(int i = 1; i < location.getCoordinates().size(); i++) {
					sb.append(",").append(location.getCoordinates().get(i).getFirst()).append(",").append(location.getCoordinates().get(i).getSecond());
				}
				
				if(id2splitKeys.containsKey(tmpReadId) && id2splitKeys.get(tmpReadId).contains(sb.toString()))
					continue;
				
				else if(writtenSplitCandidates.contains(sb.toString()))
					continue;
				
				
				writtenSplitCandidates.add(sb.toString());
				multiMappingWriter.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0\n", contextId,tmpReadId,location.getChr(),location.getStrand(),sb.toString(),location.getMismatches(),location.getStrandOfSpliceSignal(),location.overlapsKnownJunction(),mappingCounts.get(0),mappingCounts.get(1)));
				
			}
		}
	}
	
/*	private int getMismatchCount(ReadLocation location, String readSequence) {
		StringBuffer refBuffer = new StringBuffer();
		for(Pair<Integer,Integer> coordinate : location.getCoordinates()) {
			
			if(coordinate.getSecond() < coordinate.getFirst())
				return Integer.MAX_VALUE;
			
			refBuffer.append(this.currentChrSequence.substring(coordinate.getFirst() - 1, coordinate.getSecond()).toUpperCase());
		}
		
		int mismatchCount = 0;
		
		StringBuffer readBuffer = new StringBuffer();
		readBuffer.append(readSequence);
		if(location.getStrand() == '-') {
			readBuffer.reverse();
			for(int i = 0; i < readBuffer.length(); i++) {
				readBuffer.setCharAt(i,substitute(readBuffer.charAt(i)));
			}
		}
		
		
		for(int i = 0; i < readBuffer.length(); i++) {
			if(readBuffer.charAt(i) != refBuffer.charAt(i))
				mismatchCount++;
		}
		
		
		
			System.out.println(refBuffer.toString());
			System.out.println(readBuffer.toString());
			System.out.println(location.getReadId() + " "  + (mismatchCount) + " " + location.getStrand());
			System.out.println();
		
		
		return mismatchCount;
		
		
	}
*/	
	
	/**
	 * ATTENTION: The readSequence has to be in correct orientation already. 
	 * forward for alignments on the + strand and rev comp for alignments from the - strand.
	 * 
	 * 
	 * @param location
	 * @param readSequence
	 * @param localReference
	 * @param contextOffset
	 * @return
	 */
	
	private int getMismatchCountForClippedLocation(ReadLocation location, String readSequence, StringBuilder localReference, int contextOffset) {
		StringBuilder refBuffer = new StringBuilder();
		StringBuffer readBuffer = new StringBuffer();
		Pair<Integer,Integer> coordinate;
		coordinate = location.getCoordinates().get(0);
		
		//we only consider the middle segments, as the mismatch count for the first and last segment is already known.
		int readPos = coordinate.getSecond() - coordinate.getFirst() + 1;
		for(int i = 1; i < location.getCoordinates().size() - 1; i++) {
			coordinate = location.getCoordinates().get(i);
			
			if(coordinate.getFirst() > location.getCoordinates().get(i-1).getSecond()) {
				refBuffer.append(localReference.substring(coordinate.getFirst() - contextOffset, coordinate.getSecond() - contextOffset + 1).toUpperCase());
				readBuffer.append(readSequence.substring(readPos, readPos + (coordinate.getSecond() - coordinate.getFirst() + 1)));
				readPos += (coordinate.getSecond() - coordinate.getFirst() + 1);
			}
			else {
				refBuffer.append(localReference.substring(location.getCoordinates().get(i-1).getSecond() - contextOffset + 1, coordinate.getSecond() - contextOffset + 1).toUpperCase());
				readBuffer.append(readSequence.substring(readPos + location.getCoordinates().get(i-1).getSecond() - coordinate.getFirst() + 1,readPos + (coordinate.getSecond() - coordinate.getFirst() + 1)));
				readPos += (coordinate.getSecond() - coordinate.getFirst() + 1);
				
			}
		}
		
		int mismatchCount = 0;
		for(int i = 0; i < readBuffer.length(); i++) {
			if(readBuffer.charAt(i) != refBuffer.charAt(i))
				mismatchCount++;
		}
		return mismatchCount;
	}
	
	
	private int getMismatchCount(ReadLocation location, String readSequence, StringBuilder localReference, int contextOffset) {
		
		StringBuilder refBuffer = new StringBuilder();
		StringBuffer readBuffer = new StringBuffer();
		Pair<Integer,Integer> coordinate;
		
		
		if(location.getStrand() == '-') {
			readBuffer.append(readSequence);
			readBuffer.reverse();
			for(int i = 0; i < readBuffer.length(); i++) {
				readBuffer.setCharAt(i,substitute(readBuffer.charAt(i)));
			}
			readSequence = readBuffer.toString();
			readBuffer.setLength(0);
		}
		
		
		coordinate = location.getCoordinates().get(0);
		if(coordinate.getSecond() < coordinate.getFirst() || coordinate.getSecond() - contextOffset + 1 >= localReference.length()  || coordinate.getFirst() - contextOffset < 0)
			return Integer.MAX_VALUE;
		
		refBuffer.append(localReference.substring(coordinate.getFirst() - contextOffset, coordinate.getSecond() - contextOffset + 1).toUpperCase());
		readBuffer.append(readSequence.substring(readBuffer.length(), readBuffer.length() + (coordinate.getSecond() - coordinate.getFirst() + 1)));
		
		
		
		
		
		int insertionSize = 0;
		int readPos = readBuffer.length();
		for(int i = 1; i < location.getCoordinates().size(); i++) {
			coordinate = location.getCoordinates().get(i);
			if(coordinate.getSecond() < coordinate.getFirst() || coordinate.getSecond() - contextOffset + 1 >= localReference.length()  || coordinate.getFirst() - contextOffset < 0)
				return Integer.MAX_VALUE;
			

			if(coordinate.getFirst() > location.getCoordinates().get(i-1).getSecond()) {
				refBuffer.append(localReference.substring(coordinate.getFirst() - contextOffset, coordinate.getSecond() - contextOffset + 1).toUpperCase());
				readBuffer.append(readSequence.substring(readPos, readPos + (coordinate.getSecond() - coordinate.getFirst() + 1)));
				readPos += (coordinate.getSecond() - coordinate.getFirst() + 1);
			}
			else {
				refBuffer.append(localReference.substring(location.getCoordinates().get(i-1).getSecond() - contextOffset + 1, coordinate.getSecond() - contextOffset + 1).toUpperCase());
				readBuffer.append(readSequence.substring(readPos + location.getCoordinates().get(i-1).getSecond() - coordinate.getFirst() + 1,readPos + (coordinate.getSecond() - coordinate.getFirst() + 1)));
				insertionSize += location.getCoordinates().get(i-1).getSecond() - coordinate.getFirst() + 1;
				readPos += (coordinate.getSecond() - coordinate.getFirst() + 1);
				
			}
		}
		
		
		
		int mismatchCount = 0;
		
	
		
		for(int i = 0; i < readBuffer.length(); i++) {
			if(readBuffer.charAt(i) != refBuffer.charAt(i))
				mismatchCount++;
		}
		
		
		return mismatchCount;
	}

	
	private void generateNonOverlappingCoordintes(ReadLocation location) {
		ArrayList<Pair<Integer,Integer>> nonOverlappingCoordinates = new ArrayList<Pair<Integer,Integer>>();
		ArrayList<Pair<Integer,Integer>> overlappingCoordinates = location.getCoordinates();
		Pair<Integer,Integer> currentPair;
		
		nonOverlappingCoordinates.add(overlappingCoordinates.get(0));
		int start = overlappingCoordinates.get(1).getFirst();
		int stop;
		for(int i = 2; i < overlappingCoordinates.size() - 1; i+= 2) {			
			stop = overlappingCoordinates.get(i).getSecond();
			nonOverlappingCoordinates.add(new Pair<Integer,Integer>(start,stop));
			
			start = overlappingCoordinates.get(i+1).getFirst();
		}
		//add the last segment
		stop = overlappingCoordinates.get(overlappingCoordinates.size() - 1).getSecond();
		nonOverlappingCoordinates.add(new Pair<Integer,Integer>(start,stop));
		
		location.setCoordinates(nonOverlappingCoordinates);
	}
	
	
	private void generateNonOverlappingCoordintes(ReadLocation location, ArrayList<String> segmentReadSequences, String realReadSequence) {
		
		ArrayList<Pair<Integer,Integer>> nonOverlappingCoordinates = new ArrayList<Pair<Integer,Integer>>();
		ArrayList<Pair<Integer,Integer>> overlappingCoordinates = location.getCoordinates();
		ArrayList<String> tmpSegmentReadSequences = new ArrayList<String>();
		Pair<Integer,Integer> currentPair;
		StringBuffer combinedSeq = new StringBuffer();
		
		if(location.getStrand() == '+') {
			tmpSegmentReadSequences.addAll(segmentReadSequences);
		}
		else {
			for(int i = 0; i < segmentReadSequences.size(); i++) {
				combinedSeq.setLength(0);
				combinedSeq.append(segmentReadSequences.get(i));
				tmpSegmentReadSequences.add(combinedSeq.reverse().toString());
			}
		}
		
		nonOverlappingCoordinates.add(overlappingCoordinates.get(0));
		String tmpSeq = tmpSegmentReadSequences.get(0).substring(0,overlappingCoordinates.get(0).getSecond() - overlappingCoordinates.get(0).getFirst() + 1);
		
		combinedSeq.setLength(0);
		combinedSeq.append(tmpSeq);
		int start = overlappingCoordinates.get(1).getFirst();
		int overlapLength = 0;
		int stop;
		int seqIndex = 0;
		int segmentIndex = 0;
		for(int i = 2; i < overlappingCoordinates.size() - 1; i+= 2) {
			
			overlapLength = overlappingCoordinates.get(i).getFirst() - start;
			if(overlapLength >= 0) {
				tmpSeq = tmpSegmentReadSequences.get(segmentIndex).substring(tmpSeq.length(),tmpSeq.length() + overlapLength);
				combinedSeq.append(tmpSeq);
				seqIndex = 0;
			}
			else {
				seqIndex = -1 * overlapLength;
			}
			
			segmentIndex++;
			tmpSeq = tmpSegmentReadSequences.get(segmentIndex).substring(seqIndex,seqIndex + overlappingCoordinates.get(i).getSecond()-(overlappingCoordinates.get(i).getFirst() + seqIndex) + 1);
			combinedSeq.append(tmpSeq);
			
			stop = overlappingCoordinates.get(i).getSecond();
			nonOverlappingCoordinates.add(new Pair<Integer,Integer>(start,stop));
			
			start = overlappingCoordinates.get(i+1).getFirst();
		}

		//add the last segment
		tmpSeq = tmpSegmentReadSequences.get(tmpSegmentReadSequences.size()-1).substring(tmpSeq.length() + seqIndex);
		combinedSeq.append(tmpSeq);
		
		if(location.getStrand() == '-') {
			combinedSeq.reverse();
		}
		
		int offsetA = realReadSequence.indexOf(combinedSeq.toString());
		if(offsetA == - 1) {
			
			location.setCoordinates(null);
			return;
		}
		
		
		int offsetB = realReadSequence.length() - (offsetA + combinedSeq.length());
		if(location.getStrand() == '-') {
			int tmpOffset = offsetB;
			offsetB = offsetA;
			offsetA = tmpOffset;
		}
		
		stop = overlappingCoordinates.get(overlappingCoordinates.size() - 1).getSecond() + offsetB;
		nonOverlappingCoordinates.add(new Pair<Integer,Integer>(start,stop));
		Pair<Integer,Integer> modifiedStartPair = new Pair<Integer,Integer>(nonOverlappingCoordinates.get(0).getFirst()-offsetA,nonOverlappingCoordinates.get(0).getSecond());
		nonOverlappingCoordinates.set(0,modifiedStartPair);
		
		location.setCoordinates(nonOverlappingCoordinates);
	}
	
	
	

	private void determineCombinations(ReadLocation location, int currentIntervalIndex, int maxIntervalIndex, HashMap<Integer,IntervalTree<ReadLocation>> interval2tree, ArrayList<ReadLocation> combinations, int currentReadLength, String realReadSequence, String realReadSequenceRevComp, StringBuilder readSequenceBuffer, ArrayList<String> segmentReadSequences, StringBuilder localReference, int contextOffset) {
		
		IntervalTree<ReadLocation> currentTree;
		int intersectionLength;
		String intersectingReadSequence;
		ArrayList<Pair<Integer,Integer>> tmpCoordinates;
		ArrayList<String> tmpSegmentReadSequences;
		int tmpStart;
		int tmpStop;
		int tmpReadLength;
		int insertionLength;
		boolean foundValidCombination = false;
		boolean foundArtificialCombination;
		for(int i = currentIntervalIndex; i <= maxIntervalIndex; i+=2) {
			if(interval2tree.containsKey(i)) {
				currentTree = interval2tree.get(i);
				ArrayList<ReadLocation> intersectingIntervals = new ArrayList<ReadLocation>();
				intersectingIntervals = currentTree.getIntervalsIntersecting(location.getCoordinates().get(location.getCoordinates().size()-1).getFirst(), location.getCoordinates().get(location.getCoordinates().size()-1).getSecond(), intersectingIntervals);
				for(int j = 0; j < intersectingIntervals.size(); j++) {
					if(intersectingIntervals.get(j).getStrand() == location.getStrand()) {

						
						if(intersectingIntervals.get(j).getCoordinates().get(0).getFirst() < location.getCoordinates().get(location.getCoordinates().size()-2).getSecond() - this.maxMissmatches)
							continue;
						
						
						tmpStart = location.getCoordinates().get(location.getCoordinates().size()-1).getFirst();
						tmpStop = intersectingIntervals.get(j).getCoordinates().get(0).getSecond();
						
						intersectionLength = location.getCoordinates().get(location.getCoordinates().size()-1).getSecond() - location.getCoordinates().get(location.getCoordinates().size()-1).getFirst() + 1;
						tmpReadLength = currentReadLength + (tmpStop - tmpStart + 1) + (intersectingIntervals.get(j).getCoordinates().get(1).getSecond() - intersectingIntervals.get(j).getCoordinates().get(1).getFirst() + 1) - intersectionLength; 
						
						tmpSegmentReadSequences = new ArrayList<String>();
						tmpSegmentReadSequences.addAll(segmentReadSequences);
						tmpSegmentReadSequences.add(this.string2bitset.decompress(this.read2sequence.get(intersectingIntervals.get(j).getReadId())));
						
						tmpCoordinates = new ArrayList<Pair<Integer,Integer>>();
						tmpCoordinates.addAll(location.getCoordinates());
						ReadLocation tmpLocation = new ReadLocation(location.getChr(),location.getStrand(),location.getMappingType(),tmpCoordinates, location.getMismatches());
						
						tmpLocation.addCoordinate(intersectingIntervals.get(j).getCoordinates().get(0));
						tmpLocation.addCoordinate(intersectingIntervals.get(j).getCoordinates().get(1));
						tmpLocation.setStrandOfSpliceSignal(location.getStrandOfSpliceSignal());
						tmpLocation.setOverlapsKnownJunction(location.overlapsKnownJunction());
						tmpLocation.setReadId(location.getReadId());
						

						if(tmpReadLength == realReadSequence.length()) {
							
							if(tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getSecond() - contextOffset < localReference.length()) {
								
								generateNonOverlappingCoordintes(tmpLocation);
								
								
								//check if there are nested split parts, artificially created by insertions.
								foundArtificialCombination = false;
								for(int k = 1; k < tmpLocation.getCoordinates().size(); k ++) {
									if(tmpLocation.getCoordinates().get(k).getSecond() < tmpLocation.getCoordinates().get(k-1).getSecond()) {
										foundArtificialCombination = true;
										break;
									}
								}
								if(foundArtificialCombination) {
									continue;
								}
								
								int mismatches = getMismatchCount(tmpLocation,realReadSequence,localReference,contextOffset);
								
								
								if(mismatches <= this.maxMissmatches) {
									
									//check if the center blocks are >= seedLength (otherwise we found an artificial combination)
									foundArtificialCombination = false;
									for(int k = 1; k < tmpLocation.getCoordinates().size() - 1; k++) {
										if(tmpLocation.getCoordinates().get(k).getSecond() - tmpLocation.getCoordinates().get(k).getFirst() + 1 < this.seedLength) {
											foundArtificialCombination = true;
											break;
										}
									} 
									
									if(!foundArtificialCombination) {
										tmpLocation.setMismatches(mismatches);
										combinations.add(tmpLocation);
										foundValidCombination = true;
									}
								}
							}
						
						}
						
						//here we know that the current combination did not reach the actual read length yet. 
						else if(i < maxIntervalIndex && tmpReadLength < realReadSequence.length() && (tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getSecond() - tmpLocation.getCoordinates().get(0).getFirst()) < 2000000) {
							
							//check if the current combination is still valid
							foundArtificialCombination = false;
							if(tmpLocation.getCoordinates().size() >= 3) {
								int centerBlockStart = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size()-2).getFirst();
								int centerBlockStop = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getSecond();
								if(centerBlockStop - centerBlockStart + 1 < this.seedLength)
									foundArtificialCombination = true;
							}
							
							if(!foundArtificialCombination)
								determineCombinations(tmpLocation,currentIntervalIndex+2,maxIntervalIndex,interval2tree,combinations, tmpReadLength,realReadSequence,realReadSequenceRevComp,readSequenceBuffer,tmpSegmentReadSequences,localReference,contextOffset);
						}
					}
				}
			}
		}
		
		//here we have a combination spanning > 1 intron, with length < actual length
		if(! foundValidCombination && location.getCoordinates().size() > 2 && currentReadLength < realReadSequence.length() && (location.getCoordinates().get(location.getCoordinates().size() - 1).getSecond() - location.getCoordinates().get(0).getFirst()) < 2000000) {
			
			tmpCoordinates = new ArrayList<Pair<Integer,Integer>>();
			tmpCoordinates.addAll(location.getCoordinates());
			
			ReadLocation tmpLocation = new ReadLocation(location.getChr(),location.getStrand(),location.getMappingType(),tmpCoordinates, location.getMismatches());
			tmpLocation.setStrandOfSpliceSignal(location.getStrandOfSpliceSignal());
			tmpLocation.setOverlapsKnownJunction(location.overlapsKnownJunction());
			tmpLocation.setReadId(location.getReadId());
			
			generateNonOverlappingCoordintes(tmpLocation,segmentReadSequences,realReadSequence);
			if(tmpLocation.getCoordinates() != null) {
				
				//check if there are nested split parts, artificially created by insertions.
				foundArtificialCombination = false;
				for(int k = 1; k < tmpLocation.getCoordinates().size(); k ++) {
					if(tmpLocation.getCoordinates().get(k).getSecond() < tmpLocation.getCoordinates().get(k-1).getSecond()) {
						foundArtificialCombination = true;
						break;
					}
				}
				
				if(!foundArtificialCombination) {
					int mismatches = getMismatchCount(tmpLocation,realReadSequence,localReference,contextOffset);
					if(mismatches <= this.maxMissmatches) {
						//check if the center blocks are >= seedLength (otherwise we found an artificial combination)
						foundArtificialCombination = false;
						for(int k = 1; k < tmpLocation.getCoordinates().size() - 1; k++) {
							if(tmpLocation.getCoordinates().get(k).getSecond() - tmpLocation.getCoordinates().get(k).getFirst() + 1 < this.seedLength) {
								foundArtificialCombination = true;
								break;
							}
						}
						
						if(!foundArtificialCombination) {
							tmpLocation.setMismatches(mismatches);
							combinations.add(tmpLocation);
						}
					}
					
					else if(this.clipping && mismatches != Integer.MAX_VALUE) {
						int firstSegmentLength = tmpLocation.getCoordinates().get(0).getSecond() - tmpLocation.getCoordinates().get(0).getFirst() + 1; 
						ArrayList<Triplet<Integer,Integer,Integer>> clippingIndices =  getBestLocalAlignmentIndices((tmpLocation.getStrand() == '+')?realReadSequence.substring(0,firstSegmentLength):realReadSequenceRevComp.substring(0,firstSegmentLength),localReference,contextOffset, readSequenceBuffer, tmpLocation.getCoordinates().get(0).getFirst(), tmpLocation.getCoordinates().get(0).getFirst() + firstSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,true,false);
						if(!clippingIndices.isEmpty()) {
							int upstreamClippingIndex = clippingIndices.get(0).getFirst();
							int lastSegmentLength = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getSecond() - tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getFirst() + 1;
							mismatches = clippingIndices.get(0).getThird();
							clippingIndices = getBestLocalAlignmentIndices((tmpLocation.getStrand() == '+')? realReadSequence.substring(realReadSequence.length() - lastSegmentLength):realReadSequenceRevComp.substring(realReadSequence.length() - lastSegmentLength),localReference,contextOffset, readSequenceBuffer, tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getFirst(), tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size() - 1).getFirst() + lastSegmentLength -1 ,'0', 1, -4, this.maxMissmatches,1,false,true);
							if(!clippingIndices.isEmpty()) {
								
								insertionLength = 0;
								for(int i = 1; i < tmpLocation.getCoordinates().size();i++) {
									if(tmpLocation.getCoordinates().get(i).getFirst() <= tmpLocation.getCoordinates().get(i-1).getSecond()) {
										insertionLength += (-1 *(tmpLocation.getCoordinates().get(i).getFirst() - tmpLocation.getCoordinates().get(i-1).getSecond() - 1));
									}
								}
								
								int downstreamClippingIndex = realReadSequence.length() - lastSegmentLength + clippingIndices.get(0).getSecond() + insertionLength;
								downstreamClippingIndex = Math.min(realReadSequence.length() - 1, downstreamClippingIndex);
								
								mismatches += clippingIndices.get(0).getThird();
								mismatches+= getMismatchCountForClippedLocation(tmpLocation,(tmpLocation.getStrand() == '+')? realReadSequence : realReadSequenceRevComp,localReference,contextOffset);
								
								if(mismatches <= this.maxMissmatches) {
									tmpLocation.addCoordinate(new Pair<Integer,Integer>( -1 * upstreamClippingIndex, -1 * downstreamClippingIndex));
									tmpLocation.setMismatches(mismatches);
									combinations.add(tmpLocation);
									
								}
							}
						}														
					}
				}
			}
		}
	}
	                
	
	private ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>> parseAnnotatedSplitPositions(String annotationFilePath,boolean strandSpecific) {
		try {
			ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>> annotatedSpliceSites = new ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> hashedAnnotatedSplits = new HashMap<String,TreeMap<Integer,HashSet<Integer>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> hashedReverseAnnotatedSplits = new HashMap<String,TreeMap<Integer,HashSet<Integer>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentMap;
			annotatedSpliceSites.add(hashedAnnotatedSplits);
			if(strandSpecific)
				annotatedSpliceSites.add(hashedReverseAnnotatedSplits);
			
			BufferedReader br = new BufferedReader(new FileReader(new File(annotationFilePath)));
			String currentLine;
			String[] splittedLine;
			String[] splittedExonInfo;
			String[] splittedTranscriptInfo;
			String[] exonsInTranscript;
			Pattern spacePattern = Pattern.compile(" ");
			Pattern tabPattern = Pattern.compile("\t");
			Pattern commaPattern = Pattern.compile(",");
			
			String id;
			String chr;
			String strand;
			String exonInfo;
			String exonId;
			String transcriptInfo;
			String transcriptId;
			
			int start;
			int end;
			int exonStart;
			int exonEnd;
			int transcriptStart;
			int transcriptEnd;
			ArrayList<Exon> allExons = new ArrayList<Exon>();
			ArrayList<Exon> transcriptExons = new ArrayList<Exon>();
			Exon exonA;
			Exon exonB;
			
			while(br.ready()) {
				currentLine = br.readLine();
				if(currentLine.charAt(0) == '#') {
					splittedLine = spacePattern.split(currentLine);
					chr = splittedLine[2];
					strand = splittedLine[3];
					
					currentMap = hashedAnnotatedSplits;
					if(strandSpecific && strand.equals("-"))
						currentMap = hashedReverseAnnotatedSplits;
					
					
					if(!currentMap.containsKey(chr))
						currentMap.put(chr,new TreeMap<Integer,HashSet<Integer>>());
					
					currentLine = br.readLine();
					splittedLine = tabPattern.split(currentLine);
					allExons.clear();
					for(int i = 0; i < splittedLine.length; i++) {
						exonInfo = splittedLine[i];
						splittedExonInfo = spacePattern.split(exonInfo);
						exonId = splittedExonInfo[0];
						exonStart = Integer.valueOf(splittedExonInfo[1]);
						exonEnd = Integer.valueOf(splittedExonInfo[2]);
						allExons.add(new Exon(exonStart,exonEnd));
					}
					
					currentLine = br.readLine();
					splittedLine = tabPattern.split(currentLine);
					for(int i = 0; i < splittedLine.length; i++) {
						transcriptInfo = splittedLine[i];
						splittedTranscriptInfo = spacePattern.split(transcriptInfo);
						exonsInTranscript = commaPattern.split(splittedTranscriptInfo[3]);
						transcriptExons.clear();
						for(int j = 0; j < exonsInTranscript.length; j++) {
							int exonNumber = Integer.valueOf(exonsInTranscript[j]);
							transcriptExons.add(allExons.get(exonNumber));
						}
						Collections.sort(transcriptExons);
						for(int j = 0; j < transcriptExons.size() - 1; j++) {
							exonA = transcriptExons.get(j);
							exonB = transcriptExons.get(j+1);

							if(currentMap.get(chr).containsKey(exonA.getEnd()))
								currentMap.get(chr).get(exonA.getEnd()).add(exonB.getStart());
							else {
								HashSet<Integer> tmpSet = new HashSet<Integer>();
								tmpSet.add(exonB.getStart());
								currentMap.get(chr).put(exonA.getEnd(), tmpSet);
							}
							
							if(currentMap.get(chr).containsKey(exonB.getStart()))
								currentMap.get(chr).get(exonB.getStart()).add(exonA.getEnd());
							else {
								HashSet<Integer> tmpSet = new HashSet<Integer>();
								tmpSet.add(exonA.getEnd());
								currentMap.get(chr).put(exonB.getStart(), tmpSet);
							}
						}
						
					}
				}
			}
			br.close();
			return annotatedSpliceSites;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>> parseAnnotatedSplitPositionsFromGtf(String gtfFilePath, boolean strandSpecific) {
		try {
			ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>> annotatedSpliceSites = new ArrayList<HashMap<String,TreeMap<Integer,HashSet<Integer>>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> hashedAnnotatedSplits = new HashMap<String,TreeMap<Integer,HashSet<Integer>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> hashedReverseAnnotatedSplits = new HashMap<String,TreeMap<Integer,HashSet<Integer>>>();
			HashMap<String,TreeMap<Integer,HashSet<Integer>>> currentMap;
			annotatedSpliceSites.add(hashedAnnotatedSplits);
			if(strandSpecific)
				annotatedSpliceSites.add(hashedReverseAnnotatedSplits);
			
			BufferedReader br = new BufferedReader(new FileReader(new File(gtfFilePath)));
			String currentLine;
			
			while((currentLine = br.readLine()) != null) {
				if(currentLine.charAt(0) != '#')
					break;
			}
			
			//first check if we need a prefix.
			String tmpChr = this.chr2rmap.keySet().iterator().next();
			String chrPrefix = "";
			if(tmpChr.length() > 3 && tmpChr.substring(0,3).equals("chr"))
				chrPrefix = "chr";
			
			String annotationChrName = currentLine.split("\t")[0];
			if(annotationChrName.length() > 3 && annotationChrName.substring(0,3).equals("chr"))
				chrPrefix = "";
			br.close();
			
			//now go on and hash the splits.
			br = new BufferedReader(new FileReader(new File(gtfFilePath)));
			
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			Pattern semicolonPattern = Pattern.compile(";");
			String chr;
			String strand;
			int start;
			int end;
			String transcriptId;
			String prevTranscriptId = "";
			ArrayList<Exon> currentExons = new ArrayList<Exon>();
			Exon exonA;
			Exon exonB;
			while(br.ready()) {
				currentLine = br.readLine();
				
				if(currentLine.charAt(0) == '#')
					continue;
				
				splittedLine = tabPattern.split(currentLine);
				chr = chrPrefix + splittedLine[0];
				strand = splittedLine[6];
				currentMap = hashedAnnotatedSplits;
				if(strandSpecific && strand.equals("-"))
					currentMap = hashedReverseAnnotatedSplits;
				
				if(!currentMap.containsKey(chr))
					currentMap.put(chr, new TreeMap<Integer,HashSet<Integer>>());
				
				if(splittedLine[2].equals("exon")) {
					start = Integer.valueOf(splittedLine[3]);
					end = Integer.valueOf(splittedLine[4]);
					transcriptId = semicolonPattern.split(splittedLine[8])[1];
					
					if(!prevTranscriptId.equals(transcriptId)) {
						Collections.sort(currentExons);
						for(int j = 0; j < currentExons.size() - 1; j++) {
							exonA = currentExons.get(j);
							exonB = currentExons.get(j+1);

							if(currentMap.get(chr).containsKey(exonA.getEnd()))
								currentMap.get(chr).get(exonA.getEnd()).add(exonB.getStart());
							else {
								HashSet<Integer> tmpSet = new HashSet<Integer>();
								tmpSet.add(exonB.getStart());
								currentMap.get(chr).put(exonA.getEnd(), tmpSet);
							}
							
							if(currentMap.get(chr).containsKey(exonB.getStart()))
								currentMap.get(chr).get(exonB.getStart()).add(exonA.getEnd());
							else {
								HashSet<Integer> tmpSet = new HashSet<Integer>();
								tmpSet.add(exonA.getEnd());
								currentMap.get(chr).put(exonB.getStart(), tmpSet);
							}
						}
						currentExons.clear();
					}
					
					currentExons.add(new Exon(start,end));
					prevTranscriptId = transcriptId;
				}
			}
			br.close();
			return annotatedSpliceSites;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	private class Exon implements Comparable {
		private int start;
		private int end;
		
		public Exon(int start, int end) {
			this.start = start;
			this.end = end;
		}
		
		public int getStart() {
			return this.start;
		}
		
		public int getEnd() {
			return this.end;
		}

		@Override
		public int compareTo(Object o) {
			int startToCompare = ((Exon)o).start;
			return(Integer.valueOf(this.start).compareTo(startToCompare));
		}
	}
}
