package org.kidneyomics.rnaseq;

import java.util.LinkedList;
import java.util.List;

class BAMInfo {

	String sampleId;
	String bam;
	String log1;
	String log2;
	String splice1;
	String splice2;
	String dupmetrics;
	
	private BAMInfo() {
		
	}
	
	
	static List<BAMInfo> getBAMInfoFromLines(List<String> lines) {
		List<BAMInfo> res = new LinkedList<>();
		for(String line : lines) {
			BAMInfo b = getBAMInfoFromLine(line);
			res.add(b);
		}
		return res;
	}
	
	static BAMInfo getBAMInfoFromLine(String line) {
		String[] cols = line.split("\t");
		BAMInfo b = new BAMInfo();
		if(cols.length != 7) {
			throw new RuntimeException("bamlist is not properly formatted");
		} else {
			b.sampleId = cols[0];
			b.bam = cols[1];
			b.log1 = cols[2];
			b.log2 = cols[3];
			b.splice1 = cols[4];
			b.splice2 = cols[5];
			b.dupmetrics = cols[6];
		}
		return b;
	}


	public String getSampleId() {
		return sampleId;
	}


	public void setSampleId(String sampleId) {
		this.sampleId = sampleId;
	}


	public String getBam() {
		return bam;
	}


	public void setBam(String bam) {
		this.bam = bam;
	}


	public String getLog1() {
		return log1;
	}


	public void setLog1(String log1) {
		this.log1 = log1;
	}


	public String getLog2() {
		return log2;
	}


	public void setLog2(String log2) {
		this.log2 = log2;
	}


	public String getSplice1() {
		return splice1;
	}


	public void setSplice1(String splice1) {
		this.splice1 = splice1;
	}


	public String getSplice2() {
		return splice2;
	}


	public void setSplice2(String splice2) {
		this.splice2 = splice2;
	}


	public String getDupmetrics() {
		return dupmetrics;
	}


	public void setDupmetrics(String dupmetrics) {
		this.dupmetrics = dupmetrics;
	}
	
	
	
}
