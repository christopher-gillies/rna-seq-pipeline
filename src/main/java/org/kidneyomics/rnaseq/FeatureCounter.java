package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import org.kidneyomics.gtf.FeatureCount;

public interface FeatureCounter {

	void buildFeatures(File gtf, String type) throws IOException;
	
	Collection<String> getFeatureIds();
	
	FeatureCount getCounts(String featureId);
	
	List<FeatureCount> getCounts();
	
	void count(SAMRecordPair samRecordPair);
	
	double getAmbiguousReadCount();
	
	double getUnmappedReadCount();
	
	double getMappedReadCount();
	
	double getTotalCount();
	
	long getNumberOfPartiallyUnmappedReads();
	
	boolean validState();
	
	void logInfo();
}
