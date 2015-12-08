package org.kidneyomics.gtf;

import java.util.Comparator;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.kidneyomics.util.Chr2Int;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * 
 * @author cgillies
 * This implementation ignores strand
 * 
 * Returns 0 if there is overlap
 * Otherwise just uses other comparator
 * 
 * 
 * Note that this violates transitivity so it cannot be used for sorting
 * 
 * if A == B and B == C, then A == C
 * 
 * A |====|
 * 	 B  |======|
 * 		C	|====|
 * 
 * A == B and B == C but A != C
 */
public class FeatureIntersectionComparator implements Comparator<Feature> {

	private static LocationComparator comparator = new LocationComparator();
	private static Logger logger = LoggerFactory.getLogger(FeatureIntersectionComparator.class);
	
	@Override
	public int compare(Feature o1, Feature o2) {
		String myChr = o1.seqname();
		String yourChr = o2.seqname();
		
		
		int chrCmp = Integer.compare(Chr2Int.convert(myChr), Chr2Int.convert(yourChr));
		
		if(chrCmp == 0) {

			Location me = o1.location();
			Location you = o2.location();
			
			
			Location intersection = null;
			if(me.isSameStrand(you)) {
				intersection = me.intersection(you);
			} else {
				me = me.opposite();
				intersection = me.intersection(you);
			}
			
			//logger.info(intersection.bioStart() + "-" + intersection.bioEnd());
			//logger.info(intersection.length() + "");
			if(intersection != null && intersection.length() > 0) {
				return 0;
			} else {
				return comparator.compare(me, you);
			}
		} else {
			return chrCmp;
		}
	}

}
