package main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;

public class InitialRead {
	
	private String id;
	private ArrayList<InitialReadLocation> readLocations;
	private long readLocationStartPointer;
	private HashSet<String> duplicates;
	
	private int overallMappingCount;
	private int overallValidPairsCount;
	
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
	
	/**
	 * Constructor for contexts that fit into memory
	 */
	
	public InitialRead(String id, char strand, char mappingType, int startA, int endA, int missmatches, int overallMappingCount, int overallValidPairsCount) {
		this.id = id;
		this.readLocations = new ArrayList<InitialReadLocation>(2);
		
		InitialReadLocation initLocation = new InitialReadLocation(strand,mappingType,startA,endA,missmatches);
		this.readLocations.add(initLocation);
		this.duplicates = new HashSet<String>(1);
		this.overallMappingCount = overallMappingCount;
		this.overallValidPairsCount = overallValidPairsCount;
	}
	
	
	public InitialRead(String id, long pointerToFirstLocation, int overallMappingCount, int overallValidPairsCount) {
		this.id = id;
		this.readLocations = null;
		this.readLocationStartPointer = pointerToFirstLocation;
		this.duplicates = null;
		this.overallMappingCount = overallMappingCount;
		this.overallValidPairsCount = overallValidPairsCount;
		
	}
	
	public int getOverallMappingCount() {
		return this.overallMappingCount;
	}
	
	public int getOverallValidPairsCount() {
		return this.overallValidPairsCount;
	}
	
	public void addDuplicate(String id) {
		if(!this.duplicates.contains(id))
			this.duplicates.add(id);
	}
	
	
	public HashSet<String> getDuplicates() {
		return this.duplicates;
	}
	
		
	public void addLocation(char strand, char mappingType, int startA, int endA, int missmatches) {
		InitialReadLocation newLocation = new InitialReadLocation(strand,mappingType,startA,endA,missmatches);
		this.readLocations.add(newLocation);
	}
	
	public void addLocation(InitialReadLocation location) {
		this.readLocations.add(location);
	}
	
	
	public void setLocations(ArrayList<InitialReadLocation> locations) {
		this.readLocations = locations;
	}
	
	public String getId() {
		return this.id;
	}
	
	public void setId(String id) {
		this.id = id;
	}
	
	public ArrayList<InitialReadLocation> getLocations() {
		return this.readLocations;
	}
	
	public long getLocationStartPointer() {
		return this.readLocationStartPointer;
	}
	
	
	public  void removeLocation(String contextKey, InitialReadLocation location) {
		this.readLocations.remove(location);
	}


}
