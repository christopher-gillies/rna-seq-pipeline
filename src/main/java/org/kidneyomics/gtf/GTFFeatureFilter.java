package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public interface GTFFeatureFilter {
	/**
	 * 
	 * @param feature
	 * @return true if this feature passes the filter
	 */
	boolean filter(Feature feature);
}
