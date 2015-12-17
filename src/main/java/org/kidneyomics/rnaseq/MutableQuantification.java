package org.kidneyomics.rnaseq;

import org.biojava.nbio.genome.parsers.gff.Feature;

public interface MutableQuantification extends Quantification {

	void putSampleExpression(String sample, double expression);
	
	void putSampleCount(String sample, Feature expressionFeature);
	
	void putSampleRPKM(String sample, Feature expressionFeature);
	
}
