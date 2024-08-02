package main;

import java.util.concurrent.atomic.AtomicInteger;


public class IOCoordinator {

	private AtomicInteger currentInputOutputOperations;
	private AtomicInteger maxInputOutputOperations;
	
	
	
	public IOCoordinator(int maxSortingOperations) {
		this.currentInputOutputOperations = new AtomicInteger(0);
		this.maxInputOutputOperations = new AtomicInteger(maxSortingOperations);
	}
	
	
	public void inputOutputOperationStarted() {
		this.currentInputOutputOperations.incrementAndGet();
	}
	
	public void inputOutputOperationEnded() {
		this.currentInputOutputOperations.decrementAndGet();
	}
	
	public boolean hasFreeSlot() {
		return (this.currentInputOutputOperations.get() < this.maxInputOutputOperations.get());
	}
	
	
}
