package org.kidneyomics.rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Iterator;
import java.util.List;

import org.kidneyomics.rnaseq.stats.ReadPairStatistic;
import org.kidneyomics.rnaseq.stats.ReadPairStatisticsFactory;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ReadStatCalculator implements ApplicationCommand {
	
	private ApplicationOptions applicationOptions;
	
	private Logger logger;
	
	private List<ReadPairStatistic> stats;
	
	@Autowired
	public ReadStatCalculator(LoggerService loggerService, ApplicationOptions applicationOptions, ReadPairStatisticsFactory readPairStatisticsFactory) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
		stats = readPairStatisticsFactory.getBasicStatistics();
	}

	@Override
	public void doWork() throws Exception {
		String fileIn = applicationOptions.getFileIn();
		String outfile = applicationOptions.getFileOut();
		
		
		//Process BAM
		File fin = new File(fileIn);
		int pairsCounted = 0;
		try(BAMProcessor processor = BAMProcessorImplDict.getBAMProcessor(fin)) {
			
			SAMRecordPair pair = null;
			while( ( pair = processor.getNextReadPair()) != null) {
				if(pairsCounted % 100000 == 0) {
					logger.info("Read pairs read: " + pairsCounted);
				}
				
				//compute statistics
				for(ReadPairStatistic stat : stats) {
					stat.addReadPair(pair);
				}
				
				
				pairsCounted++;
				
				
			}
		}
		
		logger.info("Finished processing read pairs: " + pairsCounted + " pairs counted");
		
		File outfileRef = new File(outfile);
		if(outfileRef.exists()) {
			logger.info("Outfile already exists...deleting");
			outfileRef.delete();
		}
		
		logger.info("Writing log file to: " + outfile);
		try(BufferedWriter writer = Files.newBufferedWriter(Paths.get(outfile), Charset.defaultCharset(), StandardOpenOption.CREATE)) {
			
			//write header
			Iterator<ReadPairStatistic> headerIter = stats.iterator();
			while(headerIter.hasNext()) {
				ReadPairStatistic stat = headerIter.next();
				
				stat.appendHeader(writer);
				if(headerIter.hasNext()) {
					writer.append("\t");
				}
			}
			
			writer.append("\n");
			
			//write data
			//write header
			Iterator<ReadPairStatistic> dataIter = stats.iterator();
			while(dataIter.hasNext()) {
				ReadPairStatistic stat = dataIter.next();
				
				stat.appendStatistic(writer);
				if(dataIter.hasNext()) {
					writer.append("\t");
				}
			}
			
			writer.append("\n");
		}
		
		logger.info("Finished");
		
		
	}
	
	
	
}
