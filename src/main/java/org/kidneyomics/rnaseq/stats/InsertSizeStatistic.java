package org.kidneyomics.rnaseq.stats;

import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.SAMRecord;

class InsertSizeStatistic extends AbstractReadPairStatistic {

	/*
	 * Var(X) = E[ (X - E[X])^2 ] = E[ X^2 - 2X E[X] + E[X]^2 ] = E[ X^2 ] - 2E[X] E[X] + E[X]^2 = E[X^2] - E[X]^2
	 * this method suffers from underflow
	 */
	
	private final StreamingMeanAndVarianceCalculator meanAndVarCalc = new StreamingMeanAndVarianceCalculator();
	
	private static final String[] fields = { "MEAN_INSERT_SIZE", "SD_INSERT_SIZE" };
	
	@Override
	public void addReadPair(org.kidneyomics.rnaseq.SAMRecordPair pair) {
		meanAndVarCalc.add(pair.getMate1().getInferredInsertSize());
	}
	
	@Override
	protected void addRecord(SAMRecord record) {
		//do nothing
	}

	@Override
	public List<Double> getStatistic() {
		LinkedList<Double> res = new LinkedList<>();
		res.add(meanAndVarCalc.getMean());
		res.add(meanAndVarCalc.getSd());
		return res;
	}

	@Override
	protected String[] getFields() {
		// TODO Auto-generated method stub
		return fields;
	}

}