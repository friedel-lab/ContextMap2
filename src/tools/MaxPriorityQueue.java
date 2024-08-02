package tools;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class MaxPriorityQueue<E extends Comparable<E>> {

	private List<E> queue;
	private HashMap<E,Integer> indices;
	
	public MaxPriorityQueue(List<E> queue) {
		this.queue = queue;
		this.indices = new HashMap<E,Integer>();
		buildMaxHeap(this.queue);
	}
	
	private void buildMaxHeap(List<E> queue) {
		for(int i = queue.size() / 2 - 1; i >= 0; i--) {
			maxHeapify(queue,i,true);
		}
	}
	
	private void maxHeapify(List<E> queue, int index,boolean initialCall) {
		int largest = index;
		int l = getLeft(index);
		int r = getRight(index);
		if(l < this.queue.size() && isGreaterThan(l,largest))
			largest = l;
		if(r < this.queue.size() && isGreaterThan(r,largest))
			largest = r;
		
		if(largest != index) {
			Collections.swap(this.queue, index, largest);
			if(l < this.queue.size())
				this.indices.put(this.queue.get(l),l);
			if(r < this.queue.size())
				this.indices.put(this.queue.get(r),r);
			this.indices.put(this.queue.get(index),index);
			maxHeapify(this.queue,largest,initialCall);
		}
		
		if(initialCall) {
			if(l < this.queue.size())
				this.indices.put(this.queue.get(l),l);
			if(r < this.queue.size())
				this.indices.put(this.queue.get(r),r);
			
			this.indices.put(this.queue.get(index),index);
		}
	}
	
	
	public E getMaximum() {
		return this.queue.get(0);
	}
	
	public E extractMaximum() {
		if(this.queue.size() == 0) return null;
		E max = this.queue.get(0);
		Collections.swap(this.queue, 0, this.queue.size() - 1);
		this.indices.put(this.queue.get(0), 0);
		this.indices.remove(this.queue.get(this.queue.size() - 1));
		this.queue.remove(this.queue.size() - 1);
		maxHeapify(this.queue,0,false);
		return max;
	}
	
	
	public boolean increaseKey(E object) {
		if(!this.indices.containsKey(object))
			return false;
		int index = this.indices.get(object);
		int parentIndex;
		while(index > 0 && !isGreaterThan(getParent(index),index)) {
			parentIndex = getParent(index);
			Collections.swap(this.queue, parentIndex, index);
			this.indices.put(this.queue.get(parentIndex),parentIndex);
			this.indices.put(this.queue.get(index),index);
			index = parentIndex;
		}
		return true;
	}
	
	public boolean decreaseKey(E object) {
		if(!this.indices.containsKey(object))
			return false;
		int index = this.indices.get(object);
		maxHeapify(this.queue,index,false);
		return true;
	}
	
	//TODO
	public void insert(E e) {
		this.queue.add(e);
		this.indices.put(e, this.queue.size() - 1);
		increaseKey(e);
	}
		
	
	private int getParent(int index) {
		return((index - 1) / 2);
	}
	
	private int getLeft(int index) {
		return (2 * index) + 1;
	}
	
	private int getRight(int index) {
		return (2 * index) + 2;
	}
	
	private boolean isGreaterThan(int indexA, int indexB) {
		return(this.queue.get(indexA).compareTo(this.queue.get(indexB)) > 0);
	}
	
}
