package org.kidneyomics.rnaseq;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import org.kidneyomics.rnaseq.stats.ReadPairStatistic;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class BAMStatWriter {
	
	
	@Autowired
	ApplicationOptions applicationOptions;
	
	
	Logger logger;
	
	@Autowired
	BAMStatWriter(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	private List<ReadPairStatistic> stats = new LinkedList<>();
	
	public void writeStatFile() throws Exception {
		String fileOut = applicationOptions.getFileOut();
		String fileIn = applicationOptions.getFileIn();
		
		
		//Process BAM
		File fin = new File(fileIn);
		int pairsCounted = 0;
		try(BAMProcessor processor = BAMProcessorImplDict.getBAMProcessor(fin)) {
			
			SAMRecordPair pair = null;
			while( ( pair = processor.getNextReadPair()) != null) {

				pairsCounted++;
				if(pairsCounted % 100000 == 0) {
					logger.info("Read pairs counted: " + pairsCounted);
				}
				
			}
		}
		
	}
}
