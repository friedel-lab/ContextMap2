package main;

import java.util.ArrayList;

import augmentedTree.Interval;


public class ReadLocation implements Comparable, Interval {

	private String readId;
	private String chr;
	private char strand;
	private char mappingType;
	private char overlapsPolyAtail;
	private int startA;
	private int endA;
	private int startB;
	private int endB;
	private int mismatches;
	private double score;
	
	//'0' for no signal (needed for XS tag of SAM output)
	private char strandOfSpliceSignal;
	private boolean overlapsKnownJunction;
	private boolean isMSC;
	
	private ArrayList<Pair<Integer,Integer>> coordinates;
	
	private int hashCode;
	
	
	public ReadLocation(String chr, char strand,char mappingType, int startA, int endA, int startB, int endB, int mismatches) {
		this.chr = chr;
		this.mappingType = mappingType;
		this.strand = strand;
		this.startA = startA;
		this.endA = endA;
		this.startB = startB;
		this.endB = endB;
		this.mismatches = mismatches;
		this.score = 0;
		this.strandOfSpliceSignal = '0';
		this.overlapsKnownJunction = false;
		this.isMSC = false;
		
		this.coordinates = null;
		this.readId = null;
		
		this.hashCode = this.hashCode();
	}
	
	
	public ReadLocation(String chr, char strand,char mappingType,ArrayList<Pair<Integer,Integer>> coordinates, int missmatches) {
		this.chr = chr;
		this.mappingType = mappingType;
		this.strand = strand;
		this.coordinates = coordinates;
		this.mismatches = missmatches;
		this.score = 0;
		this.strandOfSpliceSignal = '0';
		this.overlapsKnownJunction = false;
		this.isMSC = false;
		this.readId = null;
		
		this.hashCode = this.hashCode();
	}
	
	
	public void updateLocation(String chr, char strand,char mappingType, int startA, int endA, int startB, int endB, int missmatches) {
		this.chr = chr;
		this.mappingType = mappingType;
		this.strand = strand;
		this.startA = startA;
		this.endA = endA;
		this.startB = startB;
		this.endB = endB;
		this.mismatches = missmatches;
		this.score = 0;
		this.strandOfSpliceSignal = '0';
		this.overlapsKnownJunction = false;
		this.isMSC = false;
		this.hashCode = this.hashCode();
	}
	
	public String getReadId() {
		return this.readId;
	}
	
	public void setReadId(String id) {
		this.readId = id;
		this.hashCode = this.hashCode();
	}
	

	public String getChr() {
		return chr;
	}
	

	public char getStrand() {
		return this.strand;
	}

	
	public char getMappingType() {
		return this.mappingType;
	}
	
	
	public ArrayList<Pair<Integer,Integer>> getCoordinates() {
		return this.coordinates;
	}
	
	public void setCoordinates(ArrayList<Pair<Integer,Integer>> coordinates) {
		this.coordinates = coordinates;
		this.hashCode = this.hashCode();
	}
	
	public void addCoordinate(Pair<Integer,Integer> coordinate) {
		this.coordinates.add(coordinate);
		this.hashCode = this.hashCode();
	}

	public int getStartA() {
		return startA;
	}


	public int getEndA() {
		return endA;
	}
	
	
	public void setStartB(int startB) {
		this.startB = startB;
		this.hashCode = this.hashCode();
	}
	
	public int getStartB() {
		return this.startB;
	}
	
	public void setEndB(int endB) {
		this.endB = endB;
		this.hashCode = this.hashCode();
	}
	
	public int getEndB() {
		return this.endB;
	}
	
	public int getMismatches() {
		return this.mismatches;
	}
	
	public void setMismatches(int mismatches) {
		this.mismatches = mismatches;
		this.hashCode = this.hashCode();
	}
	
	
	public void setScore(double score) {
		this.score = score;
		this.hashCode = this.hashCode();
	}
	
	public double getScore() {
		return this.score;
	}
	
	
	public boolean hasSpliceSignal() {
		return (this.strandOfSpliceSignal != '0');
	}

	public void setStrandOfSpliceSignal(char hasSpliceSignal) {
		this.strandOfSpliceSignal = hasSpliceSignal;
		this.hashCode = this.hashCode();
	}
	
	public char getStrandOfSpliceSignal() {
		return this.strandOfSpliceSignal;
	}
	
	public boolean overlapsKnownJunction() {
		return this.overlapsKnownJunction;
	}
	
	public void setOverlapsKnownJunction(boolean overlaps) {
		this.overlapsKnownJunction = overlaps;
		this.hashCode = this.hashCode();
	}
	
	public boolean isMSC() {
		return this.isMSC;
	}
	
	public void setMsc(boolean isMSC) {
		this.isMSC = isMSC;
	}
	
	public char getOverlapsPolyAtail() {
		return this.overlapsPolyAtail;
	}
	
	public void setOverlapsPolyAtail(char overlapsPolyAtail) {
		this.overlapsPolyAtail = overlapsPolyAtail;
	}
	
	

	public String getId() {
		return String.format("%s_%s_%s_%s_%s_%s", chr,strand,startA,endA,startB,endB);
	}
	
	public String printLocation() {
		return(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",chr,strand,startA,endA,startB,endB,mismatches));
	}
	
	public int compareTo(Object o) {
		ReadLocation readCompare = (ReadLocation)o;
		return(Double.valueOf(this.getScore()).compareTo(Double.valueOf(readCompare.getScore())));
	}


	@Override
	public int getStart() {
		// TODO Auto-generated method stub
		return this.coordinates.get(0).getFirst();
	}


	@Override
	public int getStop() {
		// TODO Auto-generated method stub
		return this.coordinates.get(0).getSecond();
	}
	
	
	public int hashCode() {
		if(this.hashCode != 0)
			return this.hashCode;
		
		
		int result = 17;
		if(readId != null)
			result = 31 * result + this.readId.hashCode();
		
		if(chr != null)
			result = 31 * result + this.chr.hashCode();
		
		result = 31 * result + (int)this.strand;
		result = 31 * result + (int)this.mappingType;
		result = 31 * result + startA;
		result = 31 * result + endA;
		result = 31 * result + startB;
		result = 31 * result + endB;
		result = 31 * result + mismatches;
		long tmpField = Double.doubleToLongBits(score);
		result = 31 * result + (int) (tmpField ^ (tmpField >>> 32));
		result = 31 * result + ((this.strandOfSpliceSignal == '+' || this.strandOfSpliceSignal == '-') ? 1 : 0);
		result = 31 * result + (overlapsKnownJunction ? 1 : 0);
		if(coordinates != null) {
			for(Pair<Integer,Integer> pair : coordinates) {
				result = 31 * result + pair.getFirst();
				result = 31 * result + pair.getSecond();
			}
		}
		return result;
	}
}
