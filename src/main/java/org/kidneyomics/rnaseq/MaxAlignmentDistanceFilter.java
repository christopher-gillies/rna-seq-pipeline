package org.kidneyomics.rnaseq;

public class MaxAlignmentDistanceFilter implements SAMRecordPairFilter {

	private final int max;
	
	public MaxAlignmentDistanceFilter(final int max) {
		this.max = max;
	}
	/**
	 * return true if the alignment distance is less than or equal to max
	 */
	@Override
	public boolean keep(SAMRecordPair pair) {
		Object ed1 = pair.getMate1().getAttribute("nM");
		if(ed1 == null) {
			 ed1 = pair.getMate1().getAttribute("NM");
		}
		if(ed1 != null) {
			int ed = (Integer) ed1;
			if(ed > max) {
				return false;
			}
		}
		if(pair.bothPairsAligned()) {
			Object ed2 = pair.getMate2().getAttribute("nM");
			if(ed2 == null) {
				 ed2 = pair.getMate1().getAttribute("NM");
			}
			if(ed2 != null) {
				int ed = (Integer) ed2;
				if(ed > max) {
					return false;
				}
			}
		}
		
		return true;
		
	}

}
