package org.kidneyomics.rnaseq.stats;

import java.util.LinkedList;
import java.util.List;

public class ReadPairStatisticsFactory {
	public List<ReadPairStatistic> getBasicStatistics() {
		
		LinkedList<ReadPairStatistic> stats = new LinkedList<>();
		
		stats.add( new BaseContentStatistic() );
		stats.add( new MeanBaseQualityStatistic());
		stats.add( new MeanBasesPerReadGreaterThanQ30() );
		stats.add( new InsertSizeStatistic() );
		
		return stats;
		
		
	}
}
