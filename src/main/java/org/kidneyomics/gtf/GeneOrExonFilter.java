package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class GeneOrExonFilter implements GTFFeatureFilter {

	private final static String exon = "exon";
	private final static String gene = "gene";
	
	@Override
	public boolean filter(Feature feature) {
		return feature.type().equals(exon) || feature.type().equals(gene);
	}
	
}
