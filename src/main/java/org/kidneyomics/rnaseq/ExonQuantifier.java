package org.kidneyomics.rnaseq;

import java.io.File;

import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ExonQuantifier {

	
	@Autowired
	ApplicationOptions applicationOptions;
	
	Logger logger;
	
	@Autowired
	ExonQuantifier(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	
	public void quantify() throws Exception {
		
		String gtfFile = applicationOptions.getGtf();
		String fileOut = applicationOptions.getFileOut();
		String fileIn = applicationOptions.getFileIn();
		
		//Read GTF file
		
		//Process BAM
		File fin = new File(fileIn);
		try(BAMProcessor processor = BAMProcessor.getBAMProcessor(fin)) {
			
			SAMRecordPair pair = null;
			int count = 0;
			while( ( pair = processor.getNextReadPair()) != null) {
				
			}
		}
		
		
	}
	
	
}
