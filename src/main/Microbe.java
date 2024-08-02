package main;

import augmentedTree.Interval;

public class Microbe implements Interval {
	
	private String id;
	private int start;
	private int stop;
	
	public Microbe(String id, int start, int end) {
		this.id = id;
		this.start = start;
		this.stop = end;
		
	}
	
	public String getId() {
		return this.id;
	}
	
	public int getStart() {
		return this.start;
	}
	
	public int getStop() {
		return this.stop;
	}
	
	public int getLength() {
		return this.stop - this.start;
	}
	
}
