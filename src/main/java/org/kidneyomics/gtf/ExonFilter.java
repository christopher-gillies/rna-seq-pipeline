package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class ExonFilter implements GTFFeatureFilter {

	@Override
	public boolean filter(Feature feature) {
		return feature.type().equals("exon");
	}

}
