package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public enum GTF_FORMAT {

	ENSEMBL,
	GENCODE;
	
	public static GTF_FORMAT getFormat(Feature feature) {
		if(feature.hasAttribute("gene_biotype")) {
			return ENSEMBL;
		} else {
			return GENCODE;
		}
	}
}
