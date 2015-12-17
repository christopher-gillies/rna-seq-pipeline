package org.kidneyomics.rnaseq;

import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class QuantificationFactory {

	MutableQuantification getQuantification(Feature feature, List<String> sampleIds, QuantificationType type) {
		switch(type) {
		case EXON_COUNTS:
			return ExonOrGeneQuantificationResult.getExonCountsQuantification(feature, sampleIds);
		case EXON_RPKM:
			return ExonOrGeneQuantificationResult.getExonRPKMQuantification(feature, sampleIds);
		case GENE_COUNTS:
			return ExonOrGeneQuantificationResult.getGeneCountsQuantification(feature, sampleIds);
		case GENE_RPKM:
			return ExonOrGeneQuantificationResult.getGeneRPKMQuantification(feature, sampleIds);
		default:
			throw new IllegalStateException("This opperation is not supported");
		}
	}
}
