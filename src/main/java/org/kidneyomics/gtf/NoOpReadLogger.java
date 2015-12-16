package org.kidneyomics.gtf;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;

import htsjdk.samtools.SAMRecord;

/**
 * 
 * @author cgillies
 * This class does nothing on purpose and is the default case
 * Other implementations will be setup for actual logging 
 */
public class NoOpReadLogger implements ReadLogger {

	@Override
	public void setFile(File file) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void logRead(String type, SAMRecord record, Map<Feature,Integer> mappedFeatures) {
		// TODO Auto-generated method stub
	}

	@Override
	public void close() {
		// TODO Auto-generated method stub
		
	}

}
