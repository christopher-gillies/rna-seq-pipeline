package org.kidneyomics.rnaseq;

import htsjdk.samtools.SAMRecord;

public class SAMRecordPair {
	private SAMRecord mate1;
	private SAMRecord mate2;
	
	public SAMRecordPair() {
	}

	public SAMRecord getMate1() {
		return mate1;
	}

	public void setMate1(SAMRecord mate1) {
		this.mate1 = mate1;
	}

	public SAMRecord getMate2() {
		return mate2;
	}

	public void setMate2(SAMRecord mate2) {
		this.mate2 = mate2;
	}
	
	
	
}
