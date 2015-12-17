package org.kidneyomics.rnaseq;

import java.io.IOException;

import org.biojava.nbio.genome.parsers.gff.Feature;

public interface Quantification extends Comparable<Quantification> {
	
	String printHeader();
	void appendTo(Appendable appendable) throws IOException;
	double getSampleExpression(String id);
	String getGeneId();
	String getGeneName();
	String getGeneType();
	String getChr();
	int getTranscriptionStartSite();
	int getStart();
	int getEnd();
	int getLength();
	void setLength(int length);
	char getStrand();
	QuantificationType getQuantificationType();
	String getId();
	Feature getFeature();
}
