package main;

import java.util.ArrayList;

public class ReadPair<A,B> extends Pair<A,B> implements Comparable<ReadPair> {

	private ArrayList<Pair<Integer,Integer>> validPairs;
	private Pair<Integer,Integer> topScoringPair;
	
	public ReadPair(String id) {
		super();
		this.validPairs = new ArrayList<Pair<Integer,Integer>>();
		this.topScoringPair = null;
	}
	
	public void addValidPair(Pair<Integer,Integer> pair) {
		this.validPairs.add(pair);
	}
	
	public ArrayList<Pair<Integer,Integer>> getValidPairs() {
		return this.validPairs;
	}
	
	public void clearValidPairs() {
		this.validPairs.clear();
	}
	
	public void setTopScoringPair(Pair<Integer,Integer> pair) {
		this.topScoringPair = pair;
	}

	public Pair<Integer,Integer> getTopScoringPair() {
		return this.topScoringPair;
	}

	@Override
	public int compareTo(ReadPair o) {
		ReadPair pairToCompare = (ReadPair)o;
		return(Double.valueOf(this.score).compareTo(pairToCompare.getScore()));
	}

}
