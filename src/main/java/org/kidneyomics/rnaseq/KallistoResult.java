package org.kidneyomics.rnaseq;

public class KallistoResult {

	private final String transcriptId;
	private final double length;
	private final double effectiveLength ;
	private final double counts;
	private final double tpm;
	
	public KallistoResult(String line) {
		String[] cols = line.split("\t");
		
		if(cols.length != 5) {
			throw new IllegalArgumentException("Line not formatted correctly");
		} else {
			transcriptId = cols[0];
			length = Double.parseDouble(cols[1]);
			effectiveLength = Double.parseDouble(cols[2]);
			counts = Double.parseDouble(cols[3]);
			tpm = Double.parseDouble(cols[4]);
		}
	
	}

	String getTranscriptId() {
		return transcriptId;
	}

	double getLength() {
		return length;
	}

	double getEffectiveLength() {
		return effectiveLength;
	}

	double getCounts() {
		return counts;
	}

	double getTpm() {
		return tpm;
	}
	
	
	
	
}
