package context;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.regex.Pattern;

import augmentedTree.Interval;
import augmentedTree.IntervalTree;
import main.Microbe;
import main.MutableDouble;
import main.MultiSparseReadLocation;
import main.Pair;
import main.SparseReadLocation;

public class CoverageCalculator extends Thread {
	
	private ArrayList<Integer> windowSizes;
	private int completeWindowIntervall;
	private int maxContextSize;
	
	
	private ArrayList<SparseReadLocation> startSortedLocations;
	private ArrayList<SparseReadLocation> endSortedLocations;
	private IntervalTree<Microbe> intervalTree;
	
	
	public CoverageCalculator(ArrayList<Integer> windowSizes) {
		this.windowSizes = windowSizes;
	}
	
	public CoverageCalculator(ArrayList<Integer> windowSizes, int completeWindowIntervall, int maxContextSize, ArrayList<SparseReadLocation> startSortedLocations, ArrayList<SparseReadLocation> endSortedLocations,IntervalTree<Microbe> intervalTree) {
		this.windowSizes = windowSizes;
		this.completeWindowIntervall = completeWindowIntervall;
		this.maxContextSize = maxContextSize;
		this.startSortedLocations = startSortedLocations;
		this.endSortedLocations = endSortedLocations;
		this.intervalTree = intervalTree;
	}
	
	public void run() {
		int currentWindowStart;
		int currentWindowEnd;
		int processedWindowSizes;
		
		double currentMaxCover;
		double tmpCoverage = Double.MIN_VALUE;
		double currentMaxCoverA = Double.MIN_VALUE;
		double currentMaxCoverB = Double.MIN_VALUE;
		TreeMap<Integer,MutableDouble> coverage = new TreeMap<Integer,MutableDouble>();
		NavigableMap<Integer,MutableDouble> subMap = new TreeMap<Integer,MutableDouble>();
		ArrayList<Double> windowCoverages = new ArrayList<Double>();
		SparseReadLocation tmpSparseLocation;

		int prevCoverageWindowEnd = -1;
		int currentCoverageWindowStart = -1;
		int currentCoverageWindowEnd = -1;
		
		
		HashMap<String,Double> window2coverage = new HashMap<String,Double>();
		String windowKey;
		StringBuilder keyBuilder = new StringBuilder();
		Microbe currentMicrobe = null;
		int microbeGenomeStart = -1;
		int microbeGenomeEnd = -1;
		boolean outOfGenome = false;
		for(int i = 0; i < this.startSortedLocations.size(); i++) {
			tmpSparseLocation = this.startSortedLocations.get(i);
			if(tmpSparseLocation.getAmbiguityFactor() == 0 || tmpSparseLocation.getAmbiguityFactor() == 1)
				continue;
			
			
			if(this.intervalTree != null && tmpSparseLocation.getCoordinates().get(0).getFirst() > microbeGenomeEnd) {
				if(!this.intervalTree.getIntervalsSpanning(tmpSparseLocation.getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).isEmpty())
					currentMicrobe = this.intervalTree.getIntervalsSpanning(tmpSparseLocation.getCoordinates().get(0).getFirst(),new ArrayList<Microbe>()).get(0);
				else
					currentMicrobe = this.intervalTree.getIntervalsSpanning(tmpSparseLocation.getCoordinates().get(0).getSecond(),new ArrayList<Microbe>()).get(0);
				
				microbeGenomeStart = currentMicrobe.getStart();
				microbeGenomeEnd = currentMicrobe.getStop();
			}
			
			if(window2coverage.size() > 100000)
				window2coverage.clear();
			
			windowCoverages.clear();
			
			
			if(tmpSparseLocation.getCoordinates().get(0).getFirst() > currentCoverageWindowEnd || tmpSparseLocation.getCoordinates().get(0).getSecond() + this.completeWindowIntervall > currentCoverageWindowEnd ||
			   (tmpSparseLocation.getCoordinates().size() > 1 && tmpSparseLocation.getCoordinates().get(1).getFirst() != -1 && tmpSparseLocation.getCoordinates().get(tmpSparseLocation.getCoordinates().size()-1).getFirst() + this.completeWindowIntervall > currentCoverageWindowEnd)) {
				prevCoverageWindowEnd = currentCoverageWindowEnd;
				currentCoverageWindowStart = tmpSparseLocation.getCoordinates().get(0).getFirst() - this.completeWindowIntervall;
				currentCoverageWindowEnd = tmpSparseLocation.getCoordinates().get(0).getFirst() + this.maxContextSize;
				
				
				if(prevCoverageWindowEnd > currentCoverageWindowStart) {
					subMap.clear();
					subMap.putAll(coverage.tailMap(currentCoverageWindowStart, true));
					coverage.clear();
					coverage.putAll(subMap);
					currentCoverageWindowStart = prevCoverageWindowEnd + 1;
				}
				else
					coverage.clear();
				
				updateCoverage(coverage,this.startSortedLocations,this.endSortedLocations,currentCoverageWindowStart,currentCoverageWindowEnd);
				
			}
			
			
			
			
			currentMaxCover = Double.MIN_VALUE;
			
			for(Pair<Integer,Integer> segment : tmpSparseLocation.getCoordinates()) {
				keyBuilder.setLength(0);
				windowKey = keyBuilder.append(segment.getFirst()).append("_").append(segment.getSecond()).toString();
				if(window2coverage.containsKey(windowKey)) {
					tmpCoverage = window2coverage.get(windowKey);
				}
				else {
					tmpCoverage = getMaximumReadCoverage(coverage,segment.getFirst(),segment.getSecond());
					window2coverage.put(windowKey, tmpCoverage);
				}
				
				if(tmpCoverage > currentMaxCover) {
					currentMaxCover = tmpCoverage;
				}
			}
			windowCoverages.add(currentMaxCover);
			
			
			//now go through the windows
			processedWindowSizes = 0;
			for(int j = 0; j < windowSizes.size(); j++) {
				//check upstream coverage
				 currentWindowStart = tmpSparseLocation.getCoordinates().get(0).getFirst() - processedWindowSizes - windowSizes.get(j);
				 currentWindowEnd = tmpSparseLocation.getCoordinates().get(0).getFirst() - 1 - processedWindowSizes;
				 
				 outOfGenome = false;
				 if(this.intervalTree != null) {
					 
					 //check upstream genome start
					 //in case the current window is out of the genome, we try to use the mirrored window from the downstream part of the read
					 //if this window is not valid, we use the previous valid window
					 if(currentWindowStart < microbeGenomeStart && currentWindowEnd <= microbeGenomeStart) {
						 //set the window positions to the mirrored version of the current window
						 if(tmpSparseLocation.getCoordinates().size() == 1 || tmpSparseLocation.getCoordinates().get(1).getFirst() == -1) {
								currentWindowStart = tmpSparseLocation.getCoordinates().get(0).getSecond() + 1 + processedWindowSizes;
								currentWindowEnd = tmpSparseLocation.getCoordinates().get(0).getSecond() + processedWindowSizes + windowSizes.get(j);
							}
							else {
								currentWindowStart = tmpSparseLocation.getCoordinates().get(tmpSparseLocation.getCoordinates().size() - 1).getSecond() + 1 + processedWindowSizes;
								currentWindowEnd = tmpSparseLocation.getCoordinates().get(tmpSparseLocation.getCoordinates().size() - 1).getSecond() + processedWindowSizes + windowSizes.get(j);
							}
						 
						 //mirrored window not valid use previous window coverage
						 if(currentWindowEnd > microbeGenomeEnd && currentWindowStart >= microbeGenomeEnd) {
							if(windowCoverages.size() == 1)
								currentMaxCoverA = windowCoverages.get(0);
							else
								currentMaxCoverA = windowCoverages.get(windowCoverages.size() - 2);
							
							outOfGenome = true;
						 }
						 //mirrored window overlaps with genome end
						 else if(currentWindowEnd > microbeGenomeEnd) {
							 currentWindowEnd = microbeGenomeEnd;
						 }
					 }
					 
					 //in case the window overlaps with the genome start we set the window start to genome start...
					 else if(currentWindowStart < microbeGenomeStart) {
						 currentWindowStart = microbeGenomeStart;
					 }
				 }
				 
				 if(!outOfGenome) {
					 keyBuilder.setLength(0);
					 windowKey = keyBuilder.append(currentWindowStart).append("_").append(currentWindowEnd).toString();
					 if(window2coverage.containsKey(windowKey)) {
						 currentMaxCoverA = window2coverage.get(windowKey);
					 }
					 else {
						 currentMaxCoverA = getMaximumReadCoverage(coverage,currentWindowStart,currentWindowEnd);
						 window2coverage.put(windowKey, currentMaxCoverA);
					 }
				 }
				 
				 
				 
				//check downstream coverage
				if(tmpSparseLocation.getCoordinates().size() == 1 || tmpSparseLocation.getCoordinates().get(1).getFirst() == -1) {
					currentWindowStart = tmpSparseLocation.getCoordinates().get(0).getSecond() + 1 + processedWindowSizes;
					currentWindowEnd = tmpSparseLocation.getCoordinates().get(0).getSecond() + processedWindowSizes + windowSizes.get(j);
				}
				else {
					currentWindowStart = tmpSparseLocation.getCoordinates().get(tmpSparseLocation.getCoordinates().size()-1).getSecond() + 1 + processedWindowSizes;
					currentWindowEnd = tmpSparseLocation.getCoordinates().get(tmpSparseLocation.getCoordinates().size()-1).getSecond() + processedWindowSizes + windowSizes.get(j);
				}
				
				outOfGenome = false;
				if(this.intervalTree != null) {
					//check downstream genome end
					if(currentWindowEnd > microbeGenomeEnd && currentWindowStart >= microbeGenomeEnd) {
						//set the window positions to the mirrored version of the current window
						currentWindowStart = tmpSparseLocation.getCoordinates().get(0).getFirst() - processedWindowSizes - windowSizes.get(j);
						currentWindowEnd = tmpSparseLocation.getCoordinates().get(0).getFirst() - 1 - processedWindowSizes;
						if(currentWindowStart < microbeGenomeStart && currentWindowEnd <= microbeGenomeStart) {
							currentMaxCoverB = windowCoverages.get(windowCoverages.size() - 1);
							outOfGenome = true;
						}
						else if(currentWindowStart < microbeGenomeStart) {
							currentWindowStart = microbeGenomeStart;
						}
					 }
					
					 else if(currentWindowEnd > microbeGenomeEnd) {
						 currentWindowEnd = microbeGenomeEnd;
					 }
				 }
				
				if(!outOfGenome) {
					keyBuilder.setLength(0);
					windowKey = keyBuilder.append(currentWindowStart).append("_").append(currentWindowEnd).toString();
					if(window2coverage.containsKey(windowKey)) {
						currentMaxCoverB = window2coverage.get(windowKey);
					}
					else {
						currentMaxCoverB = getMaximumReadCoverage(coverage,currentWindowStart,currentWindowEnd);
						window2coverage.put(windowKey, currentMaxCoverB);
					}
				}
				currentMaxCover = Math.max(currentMaxCoverA, currentMaxCoverB);
				//windowCoverages.add(currentMaxCover);
				windowCoverages.add(currentMaxCoverA);
				windowCoverages.add(currentMaxCoverB);
				processedWindowSizes += windowSizes.get(j);
			}
			
			tmpSparseLocation.setScore(calculateLocationScore(windowCoverages));
			
		}
	}
	
	private double calculateLocationScore(ArrayList<Double> coverages) {
		double score = 0.0;
		int windowCount = 1 + this.windowSizes.size();
		if(coverages.get(0) > 0)
			score += Math.pow(2, windowCount) * Math.log(coverages.get(0) + 1);
		windowCount--;
		for(int i = 1; i < coverages.size() - 1; i+=2) {
			if(coverages.get(i) > 0)
				score += Math.pow(2, windowCount) * Math.log(coverages.get(i) + 1);
			if(coverages.get(i+1) > 0)
				score += Math.pow(2, windowCount) * Math.log(coverages.get(i+1) + 1);
			windowCount--;
		}
		return score;
	}
	
	public double getScoreCutoff() {
		ArrayList<Double> coverages = new ArrayList<Double>();
		coverages.add(1.0);
		for(int i = 0; i < this.windowSizes.size(); i++) {
			coverages.add(1.0);
			coverages.add(1.0);
		}
		return calculateLocationScore(coverages);
	}

	
/*	private double calculateLocationScore(ArrayList<Double> coverages) {
		double score = 0.0;
		int windowCount = 1 + this.windowSizes.size();
		for(int i = 0; i < coverages.size(); i++) {
			if(coverages.get(i) > 0)
				score += Math.pow(2, windowCount) * Math.log(coverages.get(i) + 1);
			windowCount--;
		}
		
		return score;
	}
	
	private double getScoreCutoff() {
		ArrayList<Double> coverages = new ArrayList<Double>();
		coverages.add(1.0);
		for(int i = 0; i < this.windowSizes.size(); i++) {
			coverages.add(1.0);
		}
		return calculateLocationScore(coverages);
	}
	
	
*/
	
	
	
	
	
	private double getMaximumReadCoverage(TreeMap<Integer,MutableDouble> coverages, int windowStart, int windowEnd) {
		double maxCover = 0.0;
		//http://www-01.ibm.com/support/docview.wss?uid=swg1IZ47493 (java problem with submap iterator)
		try {
		NavigableMap<Integer,MutableDouble> subMap = coverages.subMap(windowStart,true, windowEnd,true);
		for(MutableDouble coverage : subMap.values()) {
			if(coverage.doubleValue() > maxCover)
				maxCover = coverage.doubleValue();
		}
		return maxCover;
		}
		catch(Exception e) {
			//e.printStackTrace();
			return maxCover;
		}
	}
	
	
	private void updateCoverage(TreeMap<Integer,MutableDouble> coverage,ArrayList<SparseReadLocation> locationsStartSorted,ArrayList<SparseReadLocation> locationsEndSorted,int windowStart, int windowEnd) {
		//parse locations with start points in the current interval
		ArrayList<Pair<Integer,Integer>> tmpCoordinates = new ArrayList<Pair<Integer,Integer>>();
		tmpCoordinates.add(new Pair<Integer,Integer>(windowEnd,windowEnd));
		tmpCoordinates.add(new Pair<Integer,Integer>(-1,-1));
		SparseReadLocation tmpLocation = new MultiSparseReadLocation(tmpCoordinates,'.');
		int insertionIndex = Collections.binarySearch(locationsStartSorted, tmpLocation,new SparseLocationStartComparator());
		
		if(insertionIndex < 0)
			insertionIndex = -1 * (insertionIndex+1);
		
		if(insertionIndex >= locationsStartSorted.size())
			insertionIndex--;
		
		
		tmpLocation = locationsStartSorted.get(insertionIndex);
		int currentStart = tmpLocation.getCoordinates().get(0).getFirst();
		int currentEnd = tmpLocation.getCoordinates().get(0).getSecond();
		double coverageWeight;
		
		
		while(insertionIndex >= 1 && currentStart > windowStart) {
			if(tmpLocation.getAmbiguityFactor() != 0) {
				coverageWeight = 1.0/(double)tmpLocation.getAmbiguityFactor();
				for(int i = currentStart; i <= currentEnd; i++) {
					if(coverage.containsKey(i)) {
						coverage.get(i).add(coverageWeight);
					}
					else
						coverage.put(i,new MutableDouble(coverageWeight));
				}

			}
			tmpLocation = locationsStartSorted.get(--insertionIndex);
			currentStart = tmpLocation.getCoordinates().get(0).getFirst();
			currentEnd = tmpLocation.getCoordinates().get(0).getSecond();
		}
		
		//now parse split locations with end points in the current interval
		if(locationsEndSorted.size() > 0) {
			tmpLocation = new MultiSparseReadLocation(tmpCoordinates,'.');
			insertionIndex = Collections.binarySearch(locationsEndSorted, tmpLocation,new SparseLocationEndComparator());
			
			if(insertionIndex < 0)
				insertionIndex = -1 * (insertionIndex+1);
			
			if(insertionIndex >= locationsEndSorted.size())
				insertionIndex--;
			
			tmpLocation = locationsEndSorted.get(insertionIndex);
			currentStart = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size()-1).getFirst();
			currentEnd = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size()-1).getSecond();
			
			while(insertionIndex >= 1 && currentStart > windowStart) {
				if(tmpLocation.getAmbiguityFactor() != 0) {
					//if(!(tmpLocation.getAmbiguityFactor() > 1 && tmpLocation.getCoordinates().size() > 1 && tmpLocation.getCoordinates().get(1).getFirst() != -1)) {
					coverageWeight = 1.0/(double)tmpLocation.getAmbiguityFactor();
						for(int i = currentStart; i <= currentEnd; i++) {
							if(coverage.containsKey(i))
								coverage.get(i).add(coverageWeight);
	
							else
								coverage.put(i,new MutableDouble(coverageWeight));
						}
					//}
				}
				
				tmpLocation = locationsEndSorted.get(--insertionIndex);
				currentStart = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size()-1).getFirst();
				currentEnd = tmpLocation.getCoordinates().get(tmpLocation.getCoordinates().size()-1).getSecond();		
			}
		}
	}
	
	private class SparseLocationStartComparator implements Comparator<SparseReadLocation> {
		public int compare(SparseReadLocation l1, SparseReadLocation l2) {
			int startA = l1.getCoordinates().get(0).getFirst();
			int startB = l2.getCoordinates().get(0).getFirst();
			return ((Integer)startA).compareTo(startB);
		}
		
	}
	
	private class SparseLocationEndComparator implements Comparator<SparseReadLocation> {
		public int compare(SparseReadLocation l1, SparseReadLocation l2) {
			int end1 = l1.getCoordinates().get(0).getSecond();
			if(l1.getCoordinates().size() > 1 && l1.getCoordinates().get(1).getFirst() != -1) 
				end1 = l1.getCoordinates().get(l1.getCoordinates().size()-1).getSecond();
			
			int end2 = l2.getCoordinates().get(0).getSecond();
			if(l2.getCoordinates().size() > 1 && l2.getCoordinates().get(1).getFirst() != -1)
				end2 = l2.getCoordinates().get(l2.getCoordinates().size()-1).getSecond();
			
			return ((Integer)end1).compareTo(end2);
		}
		
	}
	
	
	private Pair<HashMap<String,IntervalTree<SimpleInterval>>,HashMap<String,IntervalTree<SimpleInterval>>> parsePseudogenes(String gtfFilePath) {
		try {
			Pair<HashMap<String,IntervalTree<SimpleInterval>>,HashMap<String,IntervalTree<SimpleInterval>>> pseudogenes = new Pair<HashMap<String,IntervalTree<SimpleInterval>>,HashMap<String,IntervalTree<SimpleInterval>>>();
			HashMap<String,IntervalTree<SimpleInterval>> fwd_chr2pseudogenes = new HashMap<String,IntervalTree<SimpleInterval>>();
			HashMap<String,IntervalTree<SimpleInterval>> rev_chr2pseudogenes = new HashMap<String,IntervalTree<SimpleInterval>>();
			pseudogenes.setFirst(fwd_chr2pseudogenes);
			pseudogenes.setSecond(rev_chr2pseudogenes);
			
			BufferedReader br = new BufferedReader(new FileReader(new File(gtfFilePath)));
			String currentLine;
			String[] splittedLine;
			Pattern tabPattern = Pattern.compile("\t");
			
			while((currentLine = br.readLine()) != null) {
				splittedLine = tabPattern.split(currentLine);
				if(!splittedLine[2].contains("pseudo")) {
					continue;
				}
				
			}
			
			return pseudogenes;
		}
		catch(Exception e) {
			e.printStackTrace();
			return null;
		}
		
	}
	
	
	private class SimpleInterval implements Interval {

		private int start;
		private int stop;
		
		public SimpleInterval(int start, int stop) {
			this.start = start;
			this.stop = stop;
		}
		
		@Override
		public int getStart() {
			// TODO Auto-generated method stub
			return this.start;
		}

		@Override
		public int getStop() {
			// TODO Auto-generated method stub
			return this.stop;
		}
		
	}
}
