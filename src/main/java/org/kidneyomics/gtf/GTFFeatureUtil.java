package org.kidneyomics.gtf;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GTFFeatureUtil {

	private static Logger logger = LoggerFactory.getLogger(GTFFeatureUtil.class);
	/**
	 * 
	 * @param exons
	 * @return true if the exons are sorted by the start location else false
	 */
	public static boolean isSorted(List<Feature> exons) {
		if(exons.size() < 2) {
			return true;
		} else {
			Iterator<Feature> iter = exons.iterator();
			Feature current = iter.next();
			while(iter.hasNext()) {
				Feature next = iter.next();
				if(next.location().bioStart() < current.location().bioStart()) {
					return false;
				}
				current = next;
			}
			return true;
		}
	}
	
	/**
	 * 
	 * @param exons
	 * @return true if the exons are sorted by the start location else false
	 */
	public static boolean isSorted(Feature[] exons) {
		if(exons.length < 2) {
			return true;
		} else {
			for(int i = 0; i < exons.length - 1; i++) {
				Feature current = exons[i];
				Feature next = exons[i + 1];
				if(next.location().bioStart() < current.location().bioStart()) {
					return false;
				}
			}
			return true;
		}
	}
	
	/**
	 * 
	 * @param exons
	 * @return true if there is no overlap in exon features and false otherwise
	 */
	public static boolean hasNoOverlapIgnoreStrand(List<Feature> exons) {
		if(exons == null) {
			throw new IllegalArgumentException("exons cannot be null");
		}
		assert(isSorted(exons));
		
		if(exons.size() < 2) {
			return true;
		} else {
			Iterator<Feature> iter = exons.iterator();
			Feature current = iter.next();
			while(iter.hasNext()) {
				Feature next = iter.next();
				if(current.location().plus().overlaps(next.location().plus())) {
					logger.info(GTFFeatureRenderer.render(current));
					logger.info(GTFFeatureRenderer.render(next));
					return false;
				}
				current = next;
			}
			return true;
		}
	}
	
	/**
	 * 
	 * @param exons
	 * @return true if there is no overlap in exon features and false otherwise
	 */
	public static boolean hasNoOverlapIgnoreStrand(Feature[] exons) {
		if(exons == null) {
			throw new IllegalArgumentException("exons cannot be null");
		}
		assert(isSorted(exons));
		
		if(exons.length < 2) {
			return true;
		} else {
			for(int i = 0; i < exons.length - 1; i++) {
				Feature current = exons[i];
				Feature next = exons[i + 1];
				if(current.location().plus().overlaps(next.location().plus())) {
					logger.info(GTFFeatureRenderer.render(current));
					logger.info(GTFFeatureRenderer.render(next));
					return false;
				}
			}
			return true;
		}
	}
	
	private static FeatureComparator comparator = new FeatureComparator();
	
	/**
	 * 
	 * @param features to be sorted
	 */
	public static void sortFeatures(List<Feature> features) {
		if(features != null) {
			Collections.sort(features, comparator);
		}
	}
	
	/**
	 * 
	 * @param features features to be sorted
	 */
	public static void sortFeatures(Feature[] features) {
		if(features != null) {
			Arrays.sort(features, comparator);
		}
	}
	
	/**
	 * 
	 * @param Feature f1
	 * @param Feature f2
	 * @return true if these features overlap and false otherwise
	 */
	public static boolean overlapIgnoreStrand(Feature f1, Feature f2) {
		if(f1.seqname().equals(f2.seqname())) {
			Location l1 = f1.location().plus();
			Location l2 = f2.location().plus();
			return l1.overlaps(l2);
		} else {
			return false;
		}
	}
	
	/**
	 * 
	 * @param l1
	 * @param l2
	 * @return the union of the locations ignoring strand
	 */
	public static Location unionIgnoreStrand(Location l1, Location l2) {
		int l1Start = l1.bioStart();
		int l1End = l1.bioEnd();
		char l1Strand = l1.bioStrand();
		
		int l2Start = l2.bioStart();
		int l2End = l2.bioEnd();
		char l2Strand = l2.bioStrand();
		
		int unionStart = Math.min(l1Start, l2Start);
		int unionEnd = Math.max(l1End, l2End);
		
		return Location.fromBio(unionStart, unionEnd, l1Strand);
		
	}
	
	//public static double fractionOverlap(Feature)
	
}
