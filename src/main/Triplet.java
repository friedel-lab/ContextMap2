package main;

public class Triplet<A,B,C> {

	protected A firstItem;
	protected B secondItem;
	protected C thirdItem;
	protected double score;
	
	
	public Triplet() {
		this.firstItem = null;
		this.secondItem = null;
		this.thirdItem = null;
		this.score = Double.MIN_VALUE;
		
	}
	
	public Triplet(A first, B second, C third) {
		this.firstItem = first;
		this.secondItem = second;
		this.thirdItem = third;
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
	
	public C getThird() {
		return thirdItem;
	}
	
	public void setThird(C thirdItem) {
		this.thirdItem = thirdItem;
	}
	
	public double getScore() {
		return score;
	}


	public void setScore(double score) {
		this.score = score;
	}
	
	
}
