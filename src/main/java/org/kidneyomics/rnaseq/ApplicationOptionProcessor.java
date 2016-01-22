package org.kidneyomics.rnaseq;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.CommandLineRunner;
import org.springframework.stereotype.Component;
import org.springframework.util.StringUtils;

@Component
public class ApplicationOptionProcessor implements OptionProcessor {

	
	ApplicationOptions applicationOptions;
	
	
	Logger logger;
	
	
	@Autowired
	ApplicationOptionProcessor(ApplicationArguments args, LoggerService loggerService, ApplicationOptions applicationOptions) throws ParseException {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
		if(args == null) {
			throw new NullPointerException("ApplicationArguments args is null");
		}
		try {
			processInputs(args.getSourceArgs());
		} catch (Exception e) {
			logger.error(e.getMessage());
			System.exit(1);
		}
	}
	
	@Override
	public void processInputs(String[] args) throws ParseException {
		Options options = new Options();
		options.addOption("help",false,"Print the help message");
		options.addOption("noSharedMemory",false,"Do not used shared memory if there is an os issue");
		options.addOption("fastqList",true,"List of fastq files");
		options.addOption("bamList",true,"List of bam files");
		options.addOption("referenceSequence",true,"the reference sequence");
		options.addOption("outputDir",true,"the output directory");
		options.addOption("starAligner",true,"the path to star aligner");
		options.addOption("gtf",true,"the gtf file for annotations");
		options.addOption("picard",true,"the jar file for picard");
		options.addOption("readLength",true,"the read length");
		options.addOption("numberOfThreadsForGenomeIndex",true,"the number of threads to use for genome indexing");
		options.addOption("numberOfThreadsForAlign",true,"the number of threads to use for aligning");
		options.addOption("numberOfThreadsForFlux",true,"the number of threads to use for flux capacitor. Default is 2");
		options.addOption("fluxCapacitor",true,"the path to the fluxcapacitor program");
		options.addOption("countReadsInExons",false,"count all the reads in a bam file and output the results in a gtf file");
		options.addOption("computeBamStatistics",false,"compute various bam statistics for an input bam");
		
		options.addOption("mergeBamStatistics",false,"merge bam stats across samples");
		options.addOption("mergeDuplicateStatistics",false,"merge duplicate stats across samples");
		
		options.addOption("logReads",true,"file path for a log file. This file will contain information about non-mapped reads.");
		
		options.addOption("countReadsInAllSamples",false,"count all the reads for samples in a bam list and merge the results. must also specify a gtf annotation and output dir");
		
		
		options.addOption("mergeExonGtfs",false,"specify a gtf list and output directory to merge the results from countReadsInExons command");
		options.addOption("mergeExonStatFiles",false,"merge the read count statistic files from countReadsInExons command");
		
		options.addOption("mergeSTARLogs",false,"merge the statistic files from starAligner command");
		
		
		options.addOption("maxEditDistance",true,"The default for this is 6. This looks at the nH or NH tag for a read and will remove reads that have a value greater than the one specified here");
		
		options.addOption("fluxQuantifyMode",true,"This parameter is used for flux capacitor. The default is PAIRED. Here are other options: AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED");
		
		options.addOption("expressionMatrix",true,"the input file that specifies the location of an expression matrix");
		
		options.addOption("fileIn",true,"the input file for finding unique mapping reads");
		options.addOption("fileOut",true,"the output file after finding unique mapping reads");
		options.addOption("findUniqueMappedReads",false,"finds the reads that are uniquely mapped. This removes reads whose mapping quality is not 255.");
		
		options.addOption("outRPKM",false,"output the RPKM from flux or TPM from Kallisto. Default is to do counts");
		options.addOption("outGeneExpressionMatrix",false,"output a gene expression matrix from flux capacitor gtf list or kallisto out list");
		options.addOption("outTranscriptExpressionMatrix",false,"output a transcript expression matrix from flux capacitor gtf list or kallisto out list");
		options.addOption("outTranscriptRatioMatrix",false,"output a transcript ratio expression matrix from flux results or kallisto out list. For each sample you will have the transcript expression / total expression across the gene the transcript comes from");
		
		
		options.addOption("mapExpressionIds",false,"map the expression sample ids to new values. specify the fileIn option with a file containing a tab separated line for each sample that you want to remap: OLD_ID\tNEW_ID");
		
		options.addOption("kallistoMerge",false,"perform kallisto merging and not flux merge");
		options.addOption("kallisto",true,"the path to the program kallisto");
		options.addOption("kallistoThreads",true,"the number of threads for kallisto to use");
		options.addOption("referenceTranscriptome",true,"the reference transcriptome for kallisto");
		
		logger.info(StringUtils.arrayToCommaDelimitedString(args));
		
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse( options, args);
		
		
		if(cmd.getOptions().length == 0) {
			printHelp(options);
		}
		
		if(cmd.hasOption("help")) {
			printHelp(options);
			applicationOptions.setHelp(true);
		} else {
			applicationOptions.setHelp(false);
		}
		
		if(cmd.hasOption("kallistoMerge")) {
			applicationOptions.setKallistoMerge(true);
		} else {
			applicationOptions.setKallistoMerge(false);
		}
		
		if(cmd.hasOption("outRPKM")) {
			applicationOptions.setOutCounts(false);
		}
		
		if(cmd.hasOption("kallisto")) {
			applicationOptions.setKallisto(cmd.getOptionValue("kallisto"));
		}
		
		if(cmd.hasOption("kallistoThreads")) {
			applicationOptions.setNumThreadsKallisto(cmd.getOptionValue("kallistoThreads"));
		}
		
		if(cmd.hasOption("referenceTranscriptome")) {
			applicationOptions.setReferenceTranscriptome(cmd.getOptionValue("referenceTranscriptome"));
		}
		
		if(cmd.hasOption("mergeExonGtfs")) {
			applicationOptions.setMergeExonGtfs(true);
		} else {
			applicationOptions.setMergeExonGtfs(false);
		}
		
		
		if(cmd.hasOption("mergeBamStatistics")) {
			applicationOptions.setMergeBamStats(true);
		} else {
			applicationOptions.setMergeBamStats(false);
		}
		
		if(cmd.hasOption("mergeDuplicateStatistics")) {
			applicationOptions.setMergeDuplicateStats(true);
		} else {
			applicationOptions.setMergeDuplicateStats(false);
		}
		
		
		
		if(cmd.hasOption("mergeExonStatFiles")) {
			applicationOptions.setMergeExonStatFiles(true);
		} else {
			applicationOptions.setMergeExonStatFiles(false);
		}
		
		if(cmd.hasOption("mergeSTARLogs")) {
			applicationOptions.setMergeSTARLogs(true);
		} else {
			applicationOptions.setMergeSTARLogs(false);
		}
		
		if(cmd.hasOption("countReadsInAllSamples")) {
			applicationOptions.setCountReadsInAllSamples(true);
		} else {
			applicationOptions.setCountReadsInAllSamples(false);
		}
		
		
		if(cmd.hasOption("mapExpressionIds")) {
			applicationOptions.setMapExpressionIds(true);
		} else {
			applicationOptions.setMapExpressionIds(false);
		}
		
		if(cmd.hasOption("computeBamStatistics")) {
			applicationOptions.setBamStats(true);
		} else {
			applicationOptions.setBamStats(false);
		}
		
		
		if(cmd.hasOption("maxEditDistance")) {
			applicationOptions.setMaxEditDistance(Integer.parseInt(cmd.getOptionValue("maxEditDistance")));
		}
		
		if(cmd.hasOption("outGeneExpressionMatrix")) {
			applicationOptions.setOutGeneExpression(true);
		} else {
			applicationOptions.setOutGeneExpression(false);
		}
		
		if(cmd.hasOption("outTranscriptExpressionMatrix")) {
			applicationOptions.setOutTranscriptExpression(true);
		} else {
			applicationOptions.setOutTranscriptExpression(false);
		}
		
		if(cmd.hasOption("outTranscriptRatioMatrix")) {
			applicationOptions.setOutTranscriptRatioMatrix(true);
		} else {
			applicationOptions.setOutTranscriptRatioMatrix(false);
		}
		
		if(cmd.hasOption("noSharedMemory")) {
			applicationOptions.setNoSharedMemory(true);
		} else {
			applicationOptions.setNoSharedMemory(false);
		}
		
		if(cmd.hasOption("findUniqueMappedReads")) {
			applicationOptions.setFindUniqueMappedReads(true);
		} else {
			applicationOptions.setFindUniqueMappedReads(false);
		}
		
		if(cmd.hasOption("countReadsInExons")) {
			applicationOptions.setCountReadsInExons(true);
		} else {
			applicationOptions.setCountReadsInExons(false);
		}
		
		if(cmd.hasOption("starAligner")) {
			applicationOptions.setStar(cmd.getOptionValue("starAligner"));
		}
		
		if(cmd.hasOption("fluxQuantifyMode")) {
			applicationOptions.setFluxCapacitorQuantifyMode(cmd.getOptionValue("fluxQuantifyMode"));
		}
		
		if(cmd.hasOption("fluxCapacitor")) {
			applicationOptions.setFluxCapacitor(cmd.getOptionValue("fluxCapacitor"));
		}
		
		if(cmd.hasOption("referenceSequence")) {
			applicationOptions.setReferenceSequence(cmd.getOptionValue("referenceSequence"));
		}
		
		
		if(cmd.hasOption("outputDir")) {
			applicationOptions.setOutputDirectory(cmd.getOptionValue("outputDir"));
		}
		
		if(cmd.hasOption("fastqList")) {
			applicationOptions.setFastqFiles(cmd.getOptionValue("fastqList"));
		}
		
		if(cmd.hasOption("bamList")) {
			applicationOptions.setBamList(cmd.getOptionValue("bamList"));
		}
		
		if(cmd.hasOption("gtf")) {
			applicationOptions.setGtf(cmd.getOptionValue("gtf"));
		}
		
		if(cmd.hasOption("numberOfThreadsForGenomeIndex")) {
			applicationOptions.setNumThreadsGenomeIndex(cmd.getOptionValue("numberOfThreadsForGenomeIndex"));
		}
		
		if(cmd.hasOption("numberOfThreadsForFlux")) {
			applicationOptions.setNumThreadsFlux(cmd.getOptionValue("numberOfThreadsForFlux"));
		}
		
		if(cmd.hasOption("readLength")) {
			applicationOptions.setReadLength(Integer.parseInt(cmd.getOptionValue("readLength")));
		}
		
		if(cmd.hasOption("picard")) {
			applicationOptions.setPicard(cmd.getOptionValue("picard"));
		}
		
		if(cmd.hasOption("numberOfThreadsForAlign")) {
			applicationOptions.setNumThreadsAlign(cmd.getOptionValue("numberOfThreadsForAlign"));
		}
		
		
		if(cmd.hasOption("expressionMatrix")) {
			applicationOptions.setExpressionMatrix(cmd.getOptionValue("expressionMatrix"));
		}
		
		if(cmd.hasOption("fileIn")) {
			applicationOptions.setFileIn(cmd.getOptionValue("fileIn"));
		}
		
		if(cmd.hasOption("fileOut")) {
			applicationOptions.setFileOut(cmd.getOptionValue("fileOut"));
		}
		
		if(cmd.hasOption("logReads")) {
			applicationOptions.setReadLogFile(cmd.getOptionValue("logReads"));
		}
		
		
	}
	
	
	public void printHelp(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp( "RNA-seq pipeline", options );
		//System.exit(0);
	}


}
