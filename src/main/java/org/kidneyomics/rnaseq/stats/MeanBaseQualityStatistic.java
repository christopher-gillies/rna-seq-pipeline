package org.kidneyomics.rnaseq.stats;

import java.util.LinkedList;
import java.util.List;


import htsjdk.samtools.SAMRecord;

class MeanBaseQualityStatistic extends AbstractReadPairStatistic {

	private long totalBaseQuality = 0;
	private long count = 0;
	//private int offset = 33;
	private static final String[] fields = { "MEAN_PHRED_BASE_QUALITY" };
	
	@Override
	protected void addRecord(SAMRecord record) {
		byte[] qualities = record.getBaseQualities();
		
		//System.err.println("QUAL: " + record.getBaseQualityString());
		
		//the byte array is a PHRED byte array not a fastq byte array
		for(int i = 0; i < qualities.length; i++) {
			int quality = qualities[i];
			//char qualChar = (char) ((int) qualities[i] + 33);
			//System.err.println("Quality: " + qualChar + " --- " + quality + "---ASCII:" + (quality + 33));
			totalBaseQuality += quality;
			count++;
		}
		
		//System.err.println("Count: " + count);
	}
	
	@Override
	public List<Double> getStatistic() {
		LinkedList<Double> list = new LinkedList<>();
		list.add( totalBaseQuality / (double) count   );
		return list;
	}

	@Override
	protected String[] getFields() {
		return fields;
	}

	
}
