package org.kidneyomics.gtf;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;

import htsjdk.samtools.SAMRecord;

public interface ReadLogger {
	void setFile(File file);
	void logRead(String type, SAMRecord record, Map<Feature,Integer> mappedFeatures);
	void close();
}
