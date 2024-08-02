package main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Read implements Comparable {
	
	private String id;
	private double score;
	private ReadLocation topScoringLocation;
	private ArrayList<ReadLocation> readLocations;
	private ArrayList<String> duplicates;
	
	private int overallMappingCount;
	private int overallValidPairCount;
	
	/**
	 * There are three different main mapping types:
	 * F - read is mapped as a full read
	 * P - partially mapped read, either the start or the end has too many missmatches
	 * S - read is mapped as a split read
	 * 
	 * in case of a split read we initially parse it with its start and end coordinates.
	 * later we define every possible split position for this read (see extension step)
	 * 
	 * in case of a partial read, we will either redefine it as a split read or we will delete it
	 */
	
	public Read(String id) {
		this.id = id;
		this.score = 0.0;
		this.readLocations = new ArrayList<ReadLocation>(2);
		this.duplicates = new ArrayList<String>(2);
	}
	
	public Read(String id, String chr, char strand, char mappingType, int startA, int endA, int startB, int endB, int missmatches) {
		this.id = id;
		this.score = 0.0;
		this.readLocations = new ArrayList<ReadLocation>(2);
		ReadLocation initLocation = new ReadLocation(chr,strand,mappingType,startA,endA,startB,endB,missmatches);
		this.readLocations.add(initLocation);
		this.duplicates = new ArrayList<String>();
	}
	
	public void addDuplicate(String id) {
		this.duplicates.add(id);
	}
	
	public void clearDuplicates() {
		this.duplicates.clear();
	}
	
	public ArrayList<String> getDuplicates() {
		return this.duplicates;
	}
	
		
	public void addLocation(String chr, char strand, char mappingType, int startA, int endA, int startB, int endB, int missmatches) {
		ReadLocation newLocation = new ReadLocation(chr,strand,mappingType,startA,endA,startB,endB,missmatches);
		newLocation.setStartB(startB);
		newLocation.setEndB(endB);
		this.readLocations.add(newLocation);
	}
	
	public void addLocation(ReadLocation location) {
		this.readLocations.add(location);
	}
	
	public void setLocations(ArrayList<ReadLocation> locations) {
		this.readLocations = locations;
	}
	
	public String getId() {
		return this.id;
	}
	
	public void setId(String id) {
		this.id = id;
	}
	
	public double getScore() {
		return this.score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	
	public int getOverallMappingCount() {
		return this.overallMappingCount;
	}
	
	public void setOverallMappingCount(int count) {
		this.overallMappingCount = count;
	}
	
	public int getOverallValidPairCount() {
		return this.overallValidPairCount;
	}
	
	public void setOverallValidPairCount(int count) {
		this.overallValidPairCount = count;
	}
	
	public ReadLocation getTopScoringLocation() {
		return this.topScoringLocation;
	}
	
	public void setTopScoringLocation(ReadLocation location) {
		this.topScoringLocation = location;
	}
	
	public ArrayList<ReadLocation> getLocations() {
		return this.readLocations;
	}
	
	
	public  void removeLocation(String contextKey, ReadLocation location) {
		this.readLocations.remove(location);
	}
	
	@Override
	public int compareTo(Object o) {
		Read readToCompare = (Read)o;
		return(Double.valueOf(this.score).compareTo(readToCompare.getScore()));
	}
	
	


}
