package org.kidneyomics.gtf;

import java.util.Comparator;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;

public class FeatureIntersectionComparatorIgnoreStrandAndChromosome implements Comparator<Feature> {

	private static LocationComparator comparator = new LocationComparator();
	
	@Override
	public int compare(Feature o1, Feature o2) {
		Location me = o1.location();
		Location you = o2.location();
		
		
		Location intersection = null;
		if(me.isSameStrand(you)) {
			intersection = me.intersection(you);
		} else {
			me = me.opposite();
			intersection = me.intersection(you);
		}
		
		if(intersection != null && intersection.length() > 0) {
			return 0;
		} else {
			return comparator.compare(me, you);
		}
	}

}
