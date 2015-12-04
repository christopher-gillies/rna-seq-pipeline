package org.kidneyomics.gtf;

import java.util.Comparator;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.kidneyomics.util.Chr2Int;

public class FeatureComparator implements Comparator<Feature> {

	private static LocationComparator comparator = new LocationComparator();
	
	@Override
	public int compare(Feature o1, Feature o2) {
		String myChr = o1.seqname();
		String yourChr = o2.seqname();
		
		
		if(myChr.equals(yourChr)) {

			Location me = o1.location();
			Location you = o2.location();
			
			return comparator.compare(me, you);
		} else {
			return Integer.compare(Chr2Int.convert(myChr), Chr2Int.convert(yourChr));
		}
	}

}
