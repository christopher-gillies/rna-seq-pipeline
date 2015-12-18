package org.kidneyomics.rnaseq;

public interface BAMProcessor {
	SAMRecordPair getNextReadPair();
}