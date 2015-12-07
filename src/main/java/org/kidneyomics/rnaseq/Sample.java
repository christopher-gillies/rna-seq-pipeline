package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.stringtemplate.v4.ST;

public class Sample implements Comparable<Sample> {
	private String sampleId;
	private Collection<FASTQ> fastqFiles;
	private List<BAM> bamFiles;
	
	public Sample(String sampleId) {
		this.sampleId = sampleId;
		this.fastqFiles = new LinkedList<FASTQ>();
		this.bamFiles = new LinkedList<BAM>();
	}
	
	
	public void addFASTQ(FASTQ fastq) {
		this.fastqFiles.add(fastq);
	}
	
	
	public String getSampleId() {
		return sampleId;
	}


	public Collection<FASTQ> getFastqFiles() {
		return fastqFiles;
	}


	public List<BAM> getBamFiles() {
		return bamFiles;
	}

	
	public static Collection<Sample> getFastqFileList(File file) throws Exception {
		HashMap<String,Sample> samples = new HashMap<String,Sample>();
		if(!file.exists()) {
			throw new Exception("fastqFiles list does not exist!");
		}
		List<String> lines = FileUtils.readLines(file);
		for(String line : lines) {
			String[] cols = line.split("\t");
			if(cols.length < 2) {
				throw new Exception(line + "\n not formatted correctly");
			}
			
			String sampleId = cols[0];
			Sample sample = null;
			if(samples.containsKey(sampleId)) {
				sample = samples.get(sampleId);
			} else {
				sample = new Sample(sampleId);
				samples.put(sampleId, sample);
			}
			
			if(cols.length == 2) {
				FASTQ fastq1 = new FASTQ();
				fastq1.setFile1(cols[1]);
				
				if(fastq1.isValid()) {
					sample.addFASTQ(fastq1);
				} else {
					throw new Exception(fastq1.getFile1() + " not found!");
				}
				
			} else if(cols.length == 3) {
				
				FASTQ fastq1 = new FASTQ();
				fastq1.setFile1(cols[1]);
				fastq1.setFile2(cols[2]);
				
				if(fastq1.isValid()) {
					sample.addFASTQ(fastq1);
				} else {
					throw new Exception(fastq1.getFile1() + " or " + fastq1.getFile2() + " not found!");
				}
				
			} 
		}
		return samples.values();
	}
	
	public static List<Sample> getBamFileList(File file) throws Exception {
		//Assume one row per sample
		List<Sample> results = new LinkedList<Sample>();
		List<String> lines = FileUtils.readLines(file);
		for(String line : lines) {
			String[] cols = line.split("\t");
			if(cols.length < 2) {
				throw new Exception(line + "\n not formatted correctly");
			}
			
			String sampleId = cols[0];
			String bamfile = cols[1];
			Sample s = new Sample(sampleId);
			BAM bam = new BAM();
			s.bamFiles.add(bam.setBamFile(bamfile));
			
			results.add(s);
			
		}
		return results;
	}

	public String toString() {
		//ST template = new ST("<sample.sampleId>: FASTQs:<\n><sample.fastqFiles:{fastq | <fastq.file1>}>");
		ST template = new ST("SAMPLE: <sample.sampleId>:\nFASTQs:<sample.fastqFiles:{fastq | \t<fastq.file1>\t<fastq.file2>\n}>\nBAMs:<sample.bamFiles:{bam | \t<bam>}> ");
		template.add("sample", this);
		return template.render();
	}


	@Override
	public int compareTo(Sample o) {
		return this.getSampleId().compareTo(o.getSampleId());
	}
	
	
}
