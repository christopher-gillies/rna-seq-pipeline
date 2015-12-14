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
	 * @return the overlapping elements sorted
	 */
	private FeatureIntersectionComparator comparator = new FeatureIntersectionComparator();
	public List<Feature> findOverlappingFeatures(Feature[] sortedFeatures, Feature target) {
		
		assert(sortedFeatures != null);
		assert(GTFFeatureUtil.isSorted(sortedFeatures));
		assert(GTFFeatureUtil.hasNoOverlapIgnoreStrand(sortedFeatures));
		
		List<Feature> result = new LinkedList<>();
		if(target == null) {
			return result;
		}
		
		int index = Arrays.binarySearch(sortedFeatures, target, comparator);
		if(index >= 0) {
			
			result.add(sortedFeatures[index]);
			
			//Check overlap at intervals greater than the current
			for(int i = index + 1; i < sortedFeatures.length; i++) {
				if(GTFFeatureUtil.overlapIgnoreStrand(target, sortedFeatures[i])) {
					result.add(sortedFeatures[i]);
				} else {
					break;
				}
			}
			
			//Check overlap at intervals less than the current
			for(int i = index - 1; i > 0; i--) {
				if(GTFFeatureUtil.overlapIgnoreStrand(target, sortedFeatures[i])) {
					result.add(sortedFeatures[i]);
				} else {
					break;
				}
			}
			
			GTFFeatureUtil.sortFeatures(result);
			return result;
		} else {
			return result;
		}
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
