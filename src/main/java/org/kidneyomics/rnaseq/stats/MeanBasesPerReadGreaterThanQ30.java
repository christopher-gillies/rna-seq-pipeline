package org.kidneyomics.rnaseq.stats;

import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.SAMRecord;

class MeanBasesPerReadGreaterThanQ30 extends AbstractReadPairStatistic {

	private long numberOfBasesGreaterThanQ30 = 0;
	private long countOfMates = 0;
	private static final String[] fields = { "MEAN_BASES_PER_READ_GREATER_THAN_Q30" };
	
	@Override
	protected void addRecord(SAMRecord record) {
		countOfMates++;
		byte[] qualities = record.getBaseQualities();
		
		for(int i = 0; i < qualities.length; i++) {
			int quality = qualities[i];
			if(quality > 30) {
				numberOfBasesGreaterThanQ30++;
			}
		}
	}

	@Override
	public List<Double> getStatistic() {
		double result = 0.0;
		if(countOfMates > 0) {
			result = 2.0 * numberOfBasesGreaterThanQ30 / ( (double) countOfMates);
		}
		LinkedList<Double> list = new LinkedList<>();
		list.add(result);
		return list;
	}

	@Override
	protected String[] getFields() {
		return fields;
	}

}
