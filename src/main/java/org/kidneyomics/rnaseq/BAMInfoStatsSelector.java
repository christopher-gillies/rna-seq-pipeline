package org.kidneyomics.rnaseq;

public class BAMInfoStatsSelector implements BAMInfoSelector {

	@Override
	public String getField(BAMInfo info) {
		return info.getBamStats();
	}

}
