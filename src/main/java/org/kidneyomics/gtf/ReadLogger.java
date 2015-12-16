package org.kidneyomics.gtf;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

public interface ReadLogger {
	void setFile(File file);
	void logRead(String type, SAMRecord record);
	void close();
}
