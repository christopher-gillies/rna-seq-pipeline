package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.ReadLogger;

public interface FeatureCounter {

	void buildFeatures(File gtf, String type) throws IOException;
	
	Collection<String> getFeatureIds();
	
	FeatureCount getCounts(String featureId);
	
	List<Feature> getCounts();
	
	List<Feature> getGeneCounts();
	
	void count(SAMRecordPair samRecordPair);
	
	double getAmbiguousReadCount();
	
	double getUnmappedReadCount();
	
	double getMappedReadCount();
	
	double getTotalCount();
	
	long getNumberOfPartiallyUnmappedReads();
	
	boolean validState();
	
	void logInfo();

	void setReadLogger(ReadLogger readLogger);
	
	ReadLogger getReadLogger();
}
