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
		
		options.addOption("fluxQuantifyMode",true,"This paramter is used for flux capacitor. The default is PAIRED. Here are other options: AUTO, SINGLE, PAIRED, SINGLE_STRANDED, PAIRED_STRANDED");
		
		options.addOption("fileIn",true,"the input file for finding unique mapping reads");
		options.addOption("fileOut",true,"the output file after finding unique mapping reads");
		options.addOption("findUniqueMappedReads",false,"finds the reads that are uniquely mapped");
		
		options.addOption("outRPKM",false,"output the RPKM");
		options.addOption("outGeneExpressionMatrix",false,"output a gene expression matrix from flux capacitor gtf list");
		options.addOption("outTranscriptExpressionMatrix",false,"output a transcript expression matrix from flux capacitor gtf list");
		options.addOption("outTranscriptRatioMatrix",false,"output a transcript ratio expression matrix from flux results. For each sample you will have the transcript expression / total expression across the gene the transcript comes from");
		
		logger.info(StringUtils.arrayToCommaDelimitedString(args));
		
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse( options, args);
		
		
		if(cmd.getOptions().length == 0) {
			printHelp(options);
		}
		
		if(cmd.hasOption("help")) {
			printHelp(options);
		}
		
		if(cmd.hasOption("outRPKM")) {
			applicationOptions.setOutCounts(false);
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
		
		
		if(cmd.hasOption("fileIn")) {
			applicationOptions.setFileIn(cmd.getOptionValue("fileIn"));
		}
		
		if(cmd.hasOption("fileOut")) {
			applicationOptions.setFileOut(cmd.getOptionValue("fileOut"));
		}
		
		
	}
	
	
	public void printHelp(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp( "RNA-seq pipeline", options );
		//System.exit(0);
	}


}
