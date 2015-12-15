package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class GeneFilter implements GTFFeatureFilter {

	private final static String gene = "gene";
	
	@Override
	public boolean filter(Feature feature) {
		return feature.type().equals(gene);
	} 

}
