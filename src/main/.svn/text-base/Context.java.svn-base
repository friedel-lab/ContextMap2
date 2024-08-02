package main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.SortedMap;
import java.util.TreeMap;

public class Context {

	private String id;
	private String chr;
	private String strand;
	private int start;
	private int end;
	private int containedReads;
	boolean strandSpecific;
	
	private long pointerToLineOfFirstRead;
	private long pointerToLineOfLastRead;
	
	private ArrayList<InitialRead> reads;
	//key is the index of a read in the reads array
	private HashMap<String,Long> partialAndFullReads2filePointer;
	private TreeMap<Integer,MutableDouble> upstreamCoverage;
	private TreeMap<Integer,MutableDouble> downstreamCoverage;
	
	public Context(String chr, int start, int end, String strand) {
		this.chr = chr;
		this.strand = strand;
		this.start = start;
		this.end = end;
		this.id = String.format("%s_%s_%s_%s",this.chr,this.strand,this.start,this.end);
		
		this.containedReads = 0;
		if(strand.equals("+") || strand.equals("-"))
			this.strandSpecific = true;
		else
			this.strandSpecific = false;
		
		this.reads = new ArrayList<InitialRead>();
		this.partialAndFullReads2filePointer = new HashMap<String,Long>();
		
		this.upstreamCoverage = new TreeMap<Integer,MutableDouble>();
		this.downstreamCoverage = new TreeMap<Integer,MutableDouble>();
		
	}
	
	
	public void setUpstreamCoverage(SortedMap<Integer,MutableDouble> coverage) {
		this.upstreamCoverage.putAll(coverage);
	}
	
	public void setDownstreamCoverage(SortedMap<Integer,MutableDouble> coverage) {
		this.downstreamCoverage.putAll(coverage);
	}
	
	public TreeMap<Integer,MutableDouble> getUpstreamCoverage() {
		return this.upstreamCoverage;
	}
	
	public TreeMap<Integer,MutableDouble> getDownstreamCoverage() {
		return this.downstreamCoverage;
	}
	
	public String getId() {
		return this.id;
	}

	public int getStart() {
		return start;
	}

	
	public int getEnd() {
		return end;
	}

	public ArrayList<InitialRead> getReads() {
		return this.reads;
	}
	
	public InitialRead getRead(int index) {
		return this.reads.get(index);
	}
	
	
	
	public void addRead(InitialRead read) {
		this.reads.add(read);
	}
	
	public void resetReads() {
		this.reads.clear();
	}
	
	public HashMap<String,Long> getPartialAndFullReads2filePointer() {
		return this.partialAndFullReads2filePointer;
	}
	
	public void clearPartialAndFullReads2filePointer() {
		this.partialAndFullReads2filePointer.clear();
	}
	
	public int getContainedReads() {
		return this.containedReads;
	}
	
	public void setContainedReads(int containedReads) {
		this.containedReads = containedReads;
	}

	public String getChr() {
		return chr;
	}

	public String getStrand() {
		return this.strand;
	}
	
	public void setStrandSpecific(boolean strandSpecific) {
		this.strandSpecific = strandSpecific;
	}
	
	public boolean isStrandSpecific() {
		return this.strandSpecific;
	}

	public long getPointerToFirstRead() {
		return pointerToLineOfFirstRead;
	}

	public void setPointerToFirstRead(long lineOfFirstRead) {
		this.pointerToLineOfFirstRead = lineOfFirstRead;
	}

	public long getPointerToLastRead() {
		return pointerToLineOfLastRead;
	}

	public void setPointerToLastRead(long lineOfLastRead) {
		this.pointerToLineOfLastRead = lineOfLastRead;
	}
	
	public void clearUpAndDownstreamCoverages() {
		this.upstreamCoverage.clear();
		this.downstreamCoverage.clear();
	}
	
	
}
