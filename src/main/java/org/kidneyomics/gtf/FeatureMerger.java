package org.kidneyomics.gtf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;

public class FeatureMerger {

	public static List<Feature> mergeOverlappingFeaturesIgnoringStrand(List<Feature> features) {
		return mergeOverlappingFeatures(features, true);
	}
	
	public static List<Feature> mergeOverlappingFeatures(List<Feature> features) {
		return mergeOverlappingFeatures(features, false);
	}
	
	private static List<Feature> mergeOverlappingFeatures(List<Feature> features, boolean ignoreStrand) {
		
		if(features == null) {
			throw new IllegalArgumentException("features cannot be null");
		}
		
		/*
		 * validate all come from same chromosome and strand
		 */
		Feature first = features.get(0);
		for(Feature feature : features) {
			if(!first.seqname().equals(feature.seqname())) {
				throw new IllegalArgumentException("all features must be from the same chromosome");
			}
			
			if(!ignoreStrand) {
				if(first.location().bioStrand() != feature.location().bioStrand()) {
					throw new IllegalArgumentException("all features must be from the same strand");
				}
			}
		}
		
		List<Feature> mergedFeatures = new ArrayList<Feature>(features.size());
		
		if(features.size() == 0) {
			return mergedFeatures;
		}
		
		Collections.sort(features, new FeatureComparator());
		
		
		Iterator<Feature> iter = features.iterator();
		Feature current = iter.next();
		
		while(iter.hasNext()) {
			
			Feature next = iter.next();
			
			Location currentLocation = current.location();
			Location nextLocation = next.location();
			
			/*
			 * if we are ignoring the strand information, we must make sure features are compared on the same strand
			 */
			if(ignoreStrand) {
				if(!current.location().isSameStrand(next.location())) {
					nextLocation = nextLocation.opposite();
				}
			}
			
			Location intersection = currentLocation.intersection(nextLocation);
			if(intersection != null && intersection.length() > 0) {
				Location union = currentLocation.union(nextLocation);
				current = new Feature( current.seqname(),current.source(),current.type(),union,current.score(),current.frame(),current.attributes());
			} else {
				mergedFeatures.add(current);
				current = next;
			}
			
		}
		mergedFeatures.add(current);
		
		
		
		return mergedFeatures;
	}
}
