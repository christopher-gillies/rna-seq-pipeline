package org.kidneyomics.rnaseq;

public interface SAMRecordPairFilter {
	/**
	 * 
	 * @param pair
	 * @return true if this pair passes the filter
	 */
	boolean keep(SAMRecordPair pair);
}
