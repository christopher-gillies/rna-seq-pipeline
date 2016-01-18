package org.kidneyomics.rnaseq;

public class BAMInfoDuplicateSelector implements BAMInfoSelector {

	@Override
	public String getField(BAMInfo info) {
		return info.getDupmetrics();
	}

}
