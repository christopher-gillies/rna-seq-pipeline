package org.kidneyomics.rnaseq;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;

import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.ApplicationHome;
import org.springframework.stereotype.Component;
import org.springframework.util.StringUtils;

@Component
public class ApplicationOptions {


	private String star;
	private String fluxCapacitor;
	private String outputDirectory;
	private String fastqFiles;
	private String bamList;
	private String referenceSequence;
	private String numThreadsGenomeIndex = "1";
	private String numThreadsFlux = "2";
	private String gtf;
	private int readLength = 100;
	private String uncompressCommand = "zcat";
	private String picard;
	private String numThreadsAlign = "1";
	private String jarLocation;
	private String fileIn;
	private String fileOut;
	private boolean findUniqueMappedReads = false;
	private boolean noSharedMemory = false;
	private boolean outCounts = true;
	private boolean outGeneExpression;
	private boolean outTranscriptExpression;
	
	private String fluxCapacitorQuantifyMode = "PAIRED"; //--printParameters AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED
	
	Logger logger;
	
	
	@Autowired
	ApplicationOptions(LoggerService loggerService) throws UnsupportedEncodingException {
		this.logger = loggerService.getLogger(this);
		jarLocation =  new ApplicationHome(ApplicationOptions.class).getSource().getAbsolutePath();
	}
	
	public enum Mode {
		ALIGN,
		FLUX_CAPACITOR,
		ERROR,
		FIND_UNIQUE_MAPPED_READS,
		TRANSCRIPT_EXPRESSION_MATRIX,
		GENE_EXPRESSION_MATRIX,
		TRANSCRIPT_RATIO_MATRIX
		
	}
	
	
	
	public String getFluxCapacitorQuantifyMode() {
		return fluxCapacitorQuantifyMode;
	}

	public void setFluxCapacitorQuantifyMode(String fluxCapacitorQuantifyMode) {
		this.fluxCapacitorQuantifyMode = fluxCapacitorQuantifyMode;
	}

	public String getBamList() {
		return bamList;
	}

	public void setBamList(String bamList) {
		this.bamList = bamList;
	}

	public String getJarLocation() {
		return this.jarLocation;
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
	
	
	
	public String getNumThreadsFlux() {
		return numThreadsFlux;
	}

	public void setNumThreadsFlux(String numThreadsFlux) {
		this.numThreadsFlux = numThreadsFlux;
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
	
	
	
	
	public String getFluxCapacitor() {
		return fluxCapacitor;
	}

	public void setFluxCapacitor(String fluxCapacitor) {
		this.fluxCapacitor = fluxCapacitor;
	}

	public String getPicard() {
		return picard;
	}

	public void setPicard(String picard) {
		this.picard = picard;
	}
	
	

	public String getNumThreadsAlign() {
		return numThreadsAlign;
	}

	public void setNumThreadsAlign(String numThreadsAlign) {
		this.numThreadsAlign = numThreadsAlign;
	}

	public void setJarLocation(String jarLocation) {
		this.jarLocation = jarLocation;
	}

	
	
	public String getFileIn() {
		return fileIn;
	}

	public void setFileIn(String fileIn) {
		this.fileIn = fileIn;
	}

	public String getFileOut() {
		return fileOut;
	}

	public void setFileOut(String fileOut) {
		this.fileOut = fileOut;
	}

	
	
	public boolean isFindUniqueMappedReads() {
		return findUniqueMappedReads;
	}

	public void setFindUniqueMappedReads(boolean findUniqueMappedReads) {
		this.findUniqueMappedReads = findUniqueMappedReads;
	}
	
	

	public boolean isOutGeneExpression() {
		return outGeneExpression;
	}

	public void setOutGeneExpression(boolean outGeneExpression) {
		this.outGeneExpression = outGeneExpression;
	}

	public boolean isOutTranscriptExpression() {
		return outTranscriptExpression;
	}

	public void setOutTranscriptExpression(boolean outTranscriptExpression) {
		this.outTranscriptExpression = outTranscriptExpression;
	}

	public boolean isOutCounts() {
		return outCounts;
	}

	public void setOutCounts(boolean outCounts) {
		this.outCounts = outCounts;
	}

	public boolean isNoSharedMemory() {
		return noSharedMemory;
	}

	public void setNoSharedMemory(boolean noSharedMemory) {
		this.noSharedMemory = noSharedMemory;
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
			
			if(StringUtils.isEmpty(getPicard())) {
				logger.error("please specify a picard tools jar file");
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
			
		} else if(findUniqueMappedReads) {
			
			if(StringUtils.isEmpty(getFileIn())) {
				logger.error("please specify an input file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileOut())) {
				logger.error("please specify an output file");
				System.exit(1);
			}
			
			result = Mode.FIND_UNIQUE_MAPPED_READS;
			
		} else if(!StringUtils.isEmpty(getFluxCapacitor())) {

			if(StringUtils.isEmpty(getBamList())) {
				logger.error("bam files cannot be empty when performing flux capacitor");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getGtf())) {
				logger.error("please specify a gtf file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getOutputDirectory())) {
				logger.error("please specify an output");
				System.exit(1);
			}
			
			//AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED
			switch(fluxCapacitorQuantifyMode) {
			case "AUTO":
			case "SINGLE":
			case "PAIRED":
			case "SINGLE_STRANDED":
			case "PAIRED_STRANDED":
				break;
				default:
					logger.error("Please specify one of AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED for fluxQuantifyMode");
			}
			
			result = Mode.FLUX_CAPACITOR;
		} else if(outTranscriptExpression) {
			if(StringUtils.isEmpty(getGtf())) {
				logger.error("please specify a gtf file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileIn())) {
				logger.error("please specify an input file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileOut())) {
				logger.error("please specify an output file");
				System.exit(1);
			}
			
			result = Mode.TRANSCRIPT_EXPRESSION_MATRIX;
		} else if(outGeneExpression) {
			if(StringUtils.isEmpty(getGtf())) {
				logger.error("please specify a gtf file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileIn())) {
				logger.error("please specify an input file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileOut())) {
				logger.error("please specify an output file");
				System.exit(1);
			}
			
			result = Mode.GENE_EXPRESSION_MATRIX;
		}
		
		return result;
	}
	
	
}
