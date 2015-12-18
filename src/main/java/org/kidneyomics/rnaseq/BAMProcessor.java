package org.kidneyomics.rnaseq;

public interface BAMProcessor extends AutoCloseable {
	SAMRecordPair getNextReadPair();
}