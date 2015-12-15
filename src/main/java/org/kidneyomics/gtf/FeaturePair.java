package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public interface FeaturePair {
	Feature getFirst();
	Feature getSecond();
	
	boolean overlaps();
}
