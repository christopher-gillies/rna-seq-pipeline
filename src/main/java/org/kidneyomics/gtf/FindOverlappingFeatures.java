package org.kidneyomics.gtf;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class FindOverlappingFeatures {

	/**
	 * 
	 * @param longest
	 * @param sortedFeatures these should be sorted with FeatureComparator
	 * @param target
	 * @return
	 */
	private FeatureIntersectionComparator comparator = new FeatureIntersectionComparator();
	private FeatureComparator sortComparator = new FeatureComparator();
	public List<Feature> findOverlappingFeatures(int longest, Feature[] sortedFeatures, Feature target) {
		List<Feature> results = new LinkedList<Feature>();
		
		
		
		int index = Arrays.binarySearch(sortedFeatures, target, comparator);
		
		if(index  >= 0) {
			results.add(sortedFeatures[index]);
			
			//Scan to right for additional overlapping features
			// features are sorted by their start coordinates, so search for the first feature that doesn't overlap the target. Once that is found, there cannot be any more
			
			int tmpIndex = index + 1;
			while(tmpIndex < sortedFeatures.length) {
				int cmp = comparator.compare(target, sortedFeatures[tmpIndex]);
				if(cmp == 0) {
					results.add(sortedFeatures[tmpIndex]);
				} else {
					break;
				}
				tmpIndex++;
			}
			
			//Scan to the left, adding intersecting features until the distance of the starting position of the current feature to the feature stored at index is greater than the "longest" feature
			//so if then we know there cannot be any other feature overlapping the target
			//suppose f1 overlaps target and f1s start position is longer than "longest" before the target start position. 
			//Then f1 would have to be the longest feature, because it would be longer than the longest feature
			tmpIndex = index - 1;
			Feature current = sortedFeatures[tmpIndex];
			while(tmpIndex >= 0 && (current.location().bioStart() - current.location().bioStart()) <= longest) {
				int cmp = comparator.compare(target, sortedFeatures[tmpIndex]);
				if(cmp == 0) {
					results.add(sortedFeatures[tmpIndex]);
				}
				
				tmpIndex--;
				if(tmpIndex >= 0) {
					current = sortedFeatures[tmpIndex];
				}
			}
			
		}
		
		if(results.size() > 1) {
			Collections.sort(results,sortComparator);
		}
		
		return results;
	}
	
	
	/**
	 * 
	 * @param features
	 * @return the length of the longest feature
	 */
	public int findLongest(Feature[] features) {
		int longest = 0;
		for(Feature f : features) {
			int length = f.location().length();
			if(length > longest) {
				longest = length;
			}
		}
		return longest;
	}
	
}
