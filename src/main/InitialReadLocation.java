package main;

public class InitialReadLocation {

	private char strand;
	private char mappingType;
	private int startA;
	private int endA;
	private int missmatches;
	
	private boolean hasSpliceSignal;
	private boolean overlapsKnownJunction;
	
	
	private int hashCode;
	
	
	public InitialReadLocation(char strand,char mappingType, int startA, int endA, int missmatches) {
		this.mappingType = mappingType;
		this.strand = strand;
		this.startA = startA;
		this.endA = endA;
		this.missmatches = missmatches;
		this.hasSpliceSignal = false;
		this.overlapsKnownJunction = false;
		
		this.hashCode = this.hashCode();
	}

	
	public void updateLocation(char strand,char mappingType, int startA, int endA, int missmatches) {
		this.mappingType = mappingType;
		this.strand = strand;
		this.startA = startA;
		this.endA = endA;
		this.missmatches = missmatches;
		this.hasSpliceSignal = false;
		this.overlapsKnownJunction = false;
		
		this.hashCode = this.hashCode();
	}
	

	public char getStrand() {
		return this.strand;
	}

	
	public char getMappingType() {
		return this.mappingType;
	}

	public int getStartA() {
		return startA;
	}

	public int getEndA() {
		return endA;
	}
	
	
	public int getStartB() {
		return 0;
	}
	
	public int getEndB() {
		return 0;
	}
	
	public int getMissmatches() {
		return this.missmatches;
	}
	
	public void setMissmatches(int missmatches) {
		this.missmatches = missmatches;
		this.hashCode = this.hashCode();
	}
	
	
	public boolean hasSpliceSignal() {
		return hasSpliceSignal;
	}

	public void setHasSpliceSignal(boolean hasSpliceSignal) {
		this.hasSpliceSignal = hasSpliceSignal;
	}
	
	public boolean overlapsKnownJunction() {
		return this.overlapsKnownJunction;
	}
	
	public void setOverlapsKnownJunction(boolean overlaps) {
		this.overlapsKnownJunction = overlaps;
		this.hashCode = this.hashCode();
	}


	public int hashCode() {
		if(this.hashCode != 0)
			return this.hashCode;
		
		int result = 17;
		

		result = 31 * result + (int)this.strand;
		result = 31 * result + (int)this.mappingType;
		result = 31 * result + startA;
		result = 31 * result + endA;
		result = 31 * result + missmatches;
		result = 31 * result + (hasSpliceSignal ? 1 : 0);
		result = 31 * result + (overlapsKnownJunction ? 1 : 0);
		
		return result;
	}
}
