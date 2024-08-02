package main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class SparseRead {
	
	private int overallMappingCount;
	private MultiSparseReadLocation topScoringLocation;
	private ArrayList<MultiSparseReadLocation> readLocations;
	private ArrayList<Long> readLocationPointers;
	
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
	
	public SparseRead() {
		this.readLocations = new ArrayList<MultiSparseReadLocation>(2);
	}

	
	public void setOverallMappingCount(int count) {
		this.overallMappingCount = count;
	}
	
	public int getOverallMappingCount() {
		return this.overallMappingCount;
	}
	
	public void setLocations(ArrayList<MultiSparseReadLocation> locations) {
		this.readLocations = locations;
	}

	public void addLocation(MultiSparseReadLocation location) {
		this.readLocations.add(location);
	}
	
	
	public void setLocationPointers(ArrayList<Long> pointers) {
		this.readLocationPointers = pointers;
	}
	
	public void addLocationPointer(long pointer) {
		this.readLocationPointers.add(pointer);
	}
	
	
	public MultiSparseReadLocation getTopScoringLocation() {
		return this.topScoringLocation;
	}
	
	public void setTopScoringLocation(MultiSparseReadLocation location) {
		this.topScoringLocation = location;
	}
	
	public ArrayList<MultiSparseReadLocation> getLocations() {
		return this.readLocations;
	}
	
	public ArrayList<Long> getLocationPointers() {
		return this.readLocationPointers;
	}
	
	
	public  void removeLocation(String contextKey, ReadLocation location) {
		this.readLocations.remove(location);
	}


}
