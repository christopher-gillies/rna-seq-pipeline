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
	private final ModeCalculator modeCalculator = new ModeCalculator();
	
	private static final String[] fields = { "MEAN_INSERT_SIZE", "SD_INSERT_SIZE", "MODE" };
	
	@Override
	public void addReadPair(org.kidneyomics.rnaseq.SAMRecordPair pair) {
		
		pair.reorderMatesByCoordinate();
		
		int insertSize = pair.getMate2().getAlignmentStart() - pair.getMate1().getAlignmentEnd();
		
		modeCalculator.add(insertSize);
		meanAndVarCalc.add(insertSize);
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
		res.add(((Integer)modeCalculator.getMode()).doubleValue());
		return res;
	}

	@Override
	protected String[] getFields() {
		return fields;
	}

}
