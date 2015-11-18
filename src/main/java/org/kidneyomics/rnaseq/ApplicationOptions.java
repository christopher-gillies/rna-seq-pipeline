package org.kidneyomics.rnaseq;

import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.springframework.util.StringUtils;

@Component
public class ApplicationOptions {


	private String star;
	private String outputDirectory;
	private String fastqFiles;
	private String referenceSequence;
	private String numThreadsGenomeIndex;
	private String gtf;
	private int readLength = 100;
	private String uncompressCommand = "zcat";
	
	Logger logger;
	
	
	@Autowired
	ApplicationOptions(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	public enum Mode {
		ALIGN,
		ERROR
	}
	
	public String getStar() {
		return star;
	}
	public void setStar(String star) {
		this.star = star;
	}
	public String getOutputDirectory() {
		return outputDirectory;
	}
	public void setOutputDirectory(String output) {
		this.outputDirectory = output;
	}
	public String getFastqFiles() {
		return fastqFiles;
	}
	public void setFastqFiles(String fastqFiles) {
		this.fastqFiles = fastqFiles;
	}
	
	
	public String getReferenceSequence() {
		return referenceSequence;
	}
	
	public void setReferenceSequence(String referenceSequence) {
		this.referenceSequence = referenceSequence;
	}
	
	public String getNumThreadsGenomeIndex() {
		return numThreadsGenomeIndex;
	}
	
	public void setNumThreadsGenomeIndex(String numThreadsGenomeIndex) {
		this.numThreadsGenomeIndex = numThreadsGenomeIndex;
	}
	
	public String getGtf() {
		return gtf;
	}
	
	public void setGtf(String gtf) {
		this.gtf = gtf;
	}
	public Logger getLogger() {
		return logger;
	}
	public void setLogger(Logger logger) {
		this.logger = logger;
	}
	
	
	
	public int getSjdbOverhang() {
		return readLength - 1;
	}
	
	public int getReadLength() {
		return readLength;
	}
	
	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}
	
	
	
	
	public String getUncompressCommand() {
		return uncompressCommand;
	}
	public void setUncompressCommand(String uncompressCommand) {
		this.uncompressCommand = uncompressCommand;
	}
	
	public Mode validateOptions() {
		Mode result = Mode.ERROR;
		
		/*
		 * Align mode
		 */
		if(!StringUtils.isEmpty(getStar())) {
			
			if(StringUtils.isEmpty(getFastqFiles())) {
				logger.error("fastq files cannot be empty when performing star aligner");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getReferenceSequence())) {
				logger.error("Please specify a reference sequence");
				System.exit(1);
			}
			
			if(StringUtils.endsWithIgnoreCase(getReferenceSequence(), ".gz")) {
				logger.error("please uncompress reference sequence file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getOutputDirectory())) {
				logger.error("please specify an output");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getGtf())) {
				logger.error("please specify a gtf file");
				System.exit(1);
			}
			
			
			if(StringUtils.endsWithIgnoreCase(getGtf(), ".gz")) {
				logger.error("please uncompress gtf file");
				System.exit(1);
			}
			
			
			if(StringUtils.isEmpty(getNumThreadsGenomeIndex())) {
				setNumThreadsGenomeIndex("1");
			}
			
			
			
			
			result = Mode.ALIGN;
			
		}
		
		return result;
	}
	
	
}
