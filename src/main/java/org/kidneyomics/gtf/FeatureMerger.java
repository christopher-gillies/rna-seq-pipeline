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
	
	public static List<Feature> removeOverlappingFeaturesIgnoreStrand(List<Feature> features) {
		boolean ignoreStrand = true;
		if(features == null) {
			throw new IllegalArgumentException("features cannot be null");
		}
		
		/*
		 * validate all come from same chromosome and strand
		 */
		{
			Feature first = features.get(0);
			Feature current = first;
			for(Feature feature : features) {
				if(!first.seqname().equals(feature.seqname())) {
					throw new IllegalArgumentException("all features must be from the same chromosome");
				}
				
				if(!ignoreStrand) {
					if(first.location().bioStrand() != feature.location().bioStrand()) {
						throw new IllegalArgumentException("all features must be from the same strand");
					}
				}
				
				if(feature.location().bioStart() < current.location().bioStart()) {
					throw new IllegalArgumentException("Input features should be sorted");
				}
				current = feature;
			}
		}
		
		List<Feature> removedFeatures = new ArrayList<Feature>(features.size());
		
		if(features.size() == 0) {
			return removedFeatures;
		}
		
		

		Iterator<Feature> iter = features.iterator();
		Feature current = iter.next();
		boolean currentFeatureShouldBeRemoved = false;
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
				
				if(!iter.hasNext()) {
					//no more elements so remove both
					removedFeatures.add(current);
					removedFeatures.add(next);
				} else if(current.location().bioEnd() < next.location().bioEnd()) {
					//current feature is ends before next feature
					// so remove it and move on to the next feature
					currentFeatureShouldBeRemoved = true;
					removedFeatures.add(current);
					current = next;
				} else {
					//current will remain as the feature to use because it ends after the next feature
					// however we must mark it for removal and remove the next featue beacuse it will no longer be used
					currentFeatureShouldBeRemoved = true;
					removedFeatures.add(next);
				}
				
			} else if(currentFeatureShouldBeRemoved == true) {
				removedFeatures.add(current);
				currentFeatureShouldBeRemoved = false;
				current = next;
			} else {
				current = next;
			}
			
			
		}
		
		features.removeAll(removedFeatures);
		
		return removedFeatures;
		
		
	}
	
	private static List<Feature> mergeOverlappingFeatures(List<Feature> features, boolean ignoreStrand) {
		
		if(features == null) {
			throw new IllegalArgumentException("features cannot be null");
		}
		
		List<Feature> mergedFeatures = new ArrayList<Feature>(features.size());
		
		if(features.size() == 0) {
			return mergedFeatures;
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
				//System.err.println(GTFFeatureRenderer.render(current));
				//System.err.println(GTFFeatureRenderer.render(next));
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
