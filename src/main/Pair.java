package main;

import java.util.ArrayList;

public class Pair<A,B> {

	protected A firstItem;
	protected B secondItem;
	protected double score;
	
	
	
	public Pair() {
		this.firstItem = null;
		this.secondItem = null;
		this.score = Double.MIN_VALUE;
		
	}
	
	public Pair(A readA, B readB) {
		this.firstItem = readA;
		this.secondItem = readB;
	}
	

	public A getFirst() {
		return firstItem;
	}
	
	public void setFirst(A firstItem) {
		this.firstItem = firstItem;
	}


	public B getSecond() {
		return secondItem;
	}
	
	public void setSecond(B secondItem) {
		this.secondItem = secondItem;
	}
	
	public double getScore() {
		return score;
	}


	public void setScore(double score) {
		this.score = score;
	}
	
}
