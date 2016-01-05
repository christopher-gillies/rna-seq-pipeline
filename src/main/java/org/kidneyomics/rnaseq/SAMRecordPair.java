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
	
	public boolean bothPairsAligned() {
		return mate1 != null && mate2 != null;
	}
	
	public boolean isValidPair() {
		if(mate1.getMateAlignmentStart() == mate2.getAlignmentStart()
				&& mate2.getMateAlignmentStart() == mate1.getAlignmentStart() 
				&& mate1.getReadName().equals(mate2.getReadName())) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * 
	 * @param record
	 * @return true if record added and false otherwise
	 */
	public boolean addPair(SAMRecord record) {
		boolean result = true;
		if(mate1 == null) {
			mate1 = record;
		} else if(mate2 == null) {
			mate2 = record;
		} else {
			result = false;
		}
		return result;
	}
	
	//Swap mates if not in correct order
	public void reorderMates() {
		if(bothPairsAligned()) {
			SAMRecord tmp = null;
			if(!mate1.getFirstOfPairFlag()) {
				tmp = mate1;
				mate1 = mate2;
				mate2 = tmp;
			}
		}
	}
	
	
}
