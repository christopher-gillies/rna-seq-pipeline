package org.kidneyomics.rnaseq.stats;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.kidneyomics.rnaseq.SAMRecordPair;

public interface ReadPairStatistic {
	
	
	void addReadPair(SAMRecordPair pair);
	
	/**
	 * 
	 * @return the fields tab separated. THESE MUST BE IN THE SAME ORDER AS THE getStatistic() method
	 */
	String getHeader();
	
	/**
	 * 
	 * @return the a list of the statistics in the same order as the getHeader() methods
	 */
	List<Double> getStatistic();
	
	/**
	 * 
	 * @return the statistics as a map
	 */
	Map<String,Double> getStatisticAsMap();
	
	void appendStatistic(Appendable appendable) throws IOException;
	void appendHeader(Appendable appendable) throws IOException;
	
	
}
