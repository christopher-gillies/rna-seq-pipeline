package org.kidneyomics.rnaseq;

import java.io.File;

public class SampleKallistoResult implements Comparable<SampleKallistoResult> {
	private final String sampleId;
	private final File file;
	private final File dir;
	private final File stats;
	
	public SampleKallistoResult(String line) {
		String parts[] = line.split("\t");
		if(parts.length == 4) {
			sampleId = parts[0];
			dir = new File(parts[1]);
			file = new File(parts[2]);
			stats = new File(parts[3]);
		} else {
			throw new IllegalArgumentException("line not formatted correctly");
		}
	}
	
	public String getSampleId() {
		return sampleId;
	}


	public File getFile() {
		return file;
	}	

	File getDir() {
		return dir;
	}

	File getStats() {
		return stats;
	}

	@Override
	public int compareTo(SampleKallistoResult o) {
		return this.getSampleId().compareTo(o.getSampleId());
	}
}
