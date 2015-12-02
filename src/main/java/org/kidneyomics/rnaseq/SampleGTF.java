package org.kidneyomics.rnaseq;

import java.io.File;

public class SampleGTF {

	private String sampleId;
	private File file;
	
	
	public SampleGTF(String line) {
		String parts[] = line.split("\t");
		if(parts.length == 2) {
			sampleId = parts[0];
			file = new File(parts[1]);
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
	
	
	
}
