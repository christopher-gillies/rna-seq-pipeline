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
import org.springframework.stereotype.Component;
import org.springframework.util.StringUtils;

@Component
public class ApplicationOptionProcessor implements OptionProcessor {

	@Autowired
	ApplicationOptions applicationOptions;
	
	
	Logger logger;
	
	
	@Autowired
	ApplicationOptionProcessor(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	@Override
	public void processInputs(String[] args) throws ParseException {
		Options options = new Options();
		options.addOption("help",false,"Print the help message");
		options.addOption("fastqList",true,"List of fastq files");
		options.addOption("referenceSequence",true,"the reference sequence");
		options.addOption("outputDir",true,"the output directory");
		options.addOption("starAligner",true,"the path to star aligner");
		options.addOption("gtf",true,"the gtf file for annotations");
		options.addOption("readLength",true,"the read length");
		options.addOption("numberOfThreadsForGenomeIndex",true,"the number of threads to use for genome indexing");
		
		
		logger.info(StringUtils.arrayToCommaDelimitedString(args));
		
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = parser.parse( options, args);
		
		
		if(cmd.getOptions().length == 0) {
			printHelp(options);
		}
		
		if(cmd.hasOption("help")) {
			printHelp(options);
		}
		
		
		if(cmd.hasOption("starAligner")) {
			applicationOptions.setStar(cmd.getOptionValue("starAligner"));
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
		
		if(cmd.hasOption("gtf")) {
			applicationOptions.setGtf(cmd.getOptionValue("gtf"));
		}
		
		if(cmd.hasOption("numberOfThreadsForGenomeIndex")) {
			applicationOptions.setNumThreadsGenomeIndex(cmd.getOptionValue("numberOfThreadsForGenomeIndex"));
		}
		
		if(cmd.hasOption("readLength")) {
			applicationOptions.setReadLength(Integer.parseInt(cmd.getOptionValue("readLength")));
		}
		
		
	}
	
	
	public void printHelp(Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp( "RNA-seq pipeline", options );
		System.exit(0);
	}

}
