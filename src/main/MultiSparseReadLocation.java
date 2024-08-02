package main;

import java.util.ArrayList;

public class MultiSparseReadLocation implements SparseReadLocation {

	
	private int ambiguityFactor;
	private double score;
	private long filePointer;
	private char strand;
	
	private ArrayList<Pair<Integer,Integer>> coordinates;
	

	public MultiSparseReadLocation(ArrayList<Pair<Integer,Integer>> coordinates,char strand) {
		this.coordinates = coordinates;
		this.score = 0.0;
		this.ambiguityFactor = 1;
		this.strand = strand;
	}
	
	public ArrayList<Pair<Integer,Integer>> getCoordinates() {
		return this.coordinates;
	}
	
	public void setCoordinates(ArrayList<Pair<Integer,Integer>> coordinates) {
		this.coordinates = coordinates;
	}
	
	public void addCoordinate(Pair<Integer,Integer> coordinate) {
		this.coordinates.add(coordinate);
	}
	
	public double getScore() {
		return this.score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public long getFilePointer() {
		return this.filePointer;
	}
	
	public void setFilePointer(long filePointer) {
		this.filePointer = filePointer;
	}

	public int getAmbiguityFactor() {
		return ambiguityFactor;
	}

	public void setAmbiguityFactor(int ambigousLocations) {
		this.ambiguityFactor = ambigousLocations;
	}
	
	public void setStrand(char strand) {
		this.strand = strand;
	}
	
	public char getStrand() {
		return this.strand;
	}
	
}
