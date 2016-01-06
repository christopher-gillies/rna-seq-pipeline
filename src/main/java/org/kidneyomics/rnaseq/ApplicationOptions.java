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
	private String readLogFile = null;
	private String referenceSequence;
	private String numThreadsGenomeIndex = "1";
	private String numThreadsFlux = "2";
	private String gtf;
	private int readLength = 100;
	private String uncompressCommand = "zcat";
	private String picard;
	private String numThreadsAlign = "1";
	private String jarLocation;
	private String lastNonSampleColInExpressionMatrix = "strand";
	private String expressionMatrix;
	private String fileIn;
	private String fileOut;
	private boolean help = false;
	private boolean findUniqueMappedReads = false;
	private boolean noSharedMemory = false;
	private boolean outCounts = true;
	private boolean outGeneExpression = false;
	private boolean outTranscriptExpression = false;
	private boolean outTranscriptRatioMatrix = false;
	private boolean countReadsInExons = false;
	private boolean mergeExonGtfs = false;
	private boolean countReadsInAllSamples = false;
	private boolean mergeExonStatFiles = false;
	private boolean mergeSTARLogs = false;
	private boolean mapExpressionIds = false;
	private Mode mode = null;
	private int maxEditDistance = 6;
	
	private String fluxCapacitorQuantifyMode = "PAIRED"; //--printParameters AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED
	
	Logger logger;
	
	
	@Autowired
	ApplicationOptions(LoggerService loggerService) {
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
		TRANSCRIPT_RATIO_MATRIX,
		COUNT_READS_IN_EXONS,
		MERGE_EXON_COUNTS,
		MERGE_EXON_COUNTS_STATS,
		COUNT_READS_ALL_SAMPLES,
		MERGE_STAR_LOGS,
		MAP_EXPRESSION_IDS,
		HELP
	}
	
	
	
	
	public boolean isMergeSTARLogs() {
		return mergeSTARLogs;
	}

	public void setMergeSTARLogs(boolean mergeSTARLogs) {
		this.mergeSTARLogs = mergeSTARLogs;
	}

	public boolean isMergeExonStatFiles() {
		return mergeExonStatFiles;
	}

	public void setMergeExonStatFiles(boolean mergeExonStatFiles) {
		this.mergeExonStatFiles = mergeExonStatFiles;
	}

	public boolean isCountReadsInAllSamples() {
		return countReadsInAllSamples;
	}

	public void setCountReadsInAllSamples(boolean countReadsInAllSamples) {
		this.countReadsInAllSamples = countReadsInAllSamples;
	}

	public boolean isMergeExonGtfs() {
		return mergeExonGtfs;
	}

	public void setMergeExonGtfs(boolean mergeExonGtfs) {
		this.mergeExonGtfs = mergeExonGtfs;
	}

	public int getMaxEditDistance() {
		return maxEditDistance;
	}

	public void setMaxEditDistance(int maxEditDistance) {
		this.maxEditDistance = maxEditDistance;
	}

	public String getReadLogFile() {
		return readLogFile;
	}

	public void setReadLogFile(String readLogFile) {
		this.readLogFile = readLogFile;
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

	
	
	public boolean isOutTranscriptRatioMatrix() {
		return outTranscriptRatioMatrix;
	}

	public void setOutTranscriptRatioMatrix(boolean outTranscriptRatioMatrix) {
		this.outTranscriptRatioMatrix = outTranscriptRatioMatrix;
	}

	
	public boolean isCountReadsInExons() {
		return countReadsInExons;
	}

	public void setCountReadsInExons(boolean countReadsInExons) {
		this.countReadsInExons = countReadsInExons;
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
		} else if(outTranscriptRatioMatrix) {
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
			
			result = Mode.TRANSCRIPT_RATIO_MATRIX;
		} else if(countReadsInExons) {
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
			
			result = Mode.COUNT_READS_IN_EXONS;
		} else if(mergeExonGtfs) {
			
			if(StringUtils.isEmpty(getFileIn())) {
				logger.error("please specify an input file");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getOutputDirectory())) {
				logger.error("please specify an output directory");
				System.exit(1);
			}
			
			result = Mode.MERGE_EXON_COUNTS;
		} else if(countReadsInAllSamples) {
			
			if(StringUtils.isEmpty(getBamList())) {
				logger.error("bam files cannot be empty when performing exon counting across samples");
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
			
			result = Mode.COUNT_READS_ALL_SAMPLES;
		} else if(mergeExonStatFiles) {
			
			if(StringUtils.isEmpty(getFileIn())) {
				logger.error("please specify an input file. This should be a list of stat files. Two columns [ID]\t[STAT_FILE]");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileOut())) {
				logger.error("please specify an output file. This is the merged matrix.");
				System.exit(1);
			}
			
			result = Mode.MERGE_EXON_COUNTS_STATS;
		} else if(mergeSTARLogs) {
			
			if(StringUtils.isEmpty(getBamList())) {
				logger.error("bam files cannot be empty when performing exon counting across samples");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileOut())) {
				logger.error("please specify an output file. This is the merged matrix.");
				System.exit(1);
			}
			
			result = Mode.MERGE_STAR_LOGS;
		} else if(mapExpressionIds) {
			
			if(StringUtils.isEmpty(getFileIn())) {
				logger.error("please specify an input file. This should be a list of old its and new ids. Two columns [OLD_ID]\t[NEW_ID]");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getExpressionMatrix())) {
				logger.error("the expression matrix that you wish to remap ids for");
				System.exit(1);
			}
			
			if(StringUtils.isEmpty(getFileOut())) {
				logger.error("please specify an output file. The expression matrix with new sample ids");
				System.exit(1);
			}
			
			result = Mode.MAP_EXPRESSION_IDS;
		} else if(help) {
			result = Mode.HELP;
		}
		
		this.mode = result;
		return result;
	}

	public Mode getMode() {
		return mode;
	}

	public void setMode(Mode mode) {
		this.mode = mode;
	}

	public String getExpressionMatrix() {
		return expressionMatrix;
	}

	public void setExpressionMatrix(String expressionMatrix) {
		this.expressionMatrix = expressionMatrix;
	}

	public String getLastNonSampleColInExpressionMatrix() {
		return lastNonSampleColInExpressionMatrix;
	}

	public void setLastNonSampleColInExpressionMatrix(String lastNonSampleColInExpressionMatrix) {
		this.lastNonSampleColInExpressionMatrix = lastNonSampleColInExpressionMatrix;
	}

	public boolean isMapExpressionIds() {
		return mapExpressionIds;
	}

	public void setMapExpressionIds(boolean mapExpressionIds) {
		this.mapExpressionIds = mapExpressionIds;
	}

	public boolean isHelp() {
		return help;
	}

	public void setHelp(boolean help) {
		this.help = help;
	}
	
	
	
	
}
