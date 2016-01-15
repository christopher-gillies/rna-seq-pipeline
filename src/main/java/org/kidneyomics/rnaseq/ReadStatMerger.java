package org.kidneyomics.rnaseq;

import java.io.File;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.kidneyomics.rnaseq.stats.ReadPairStatisticsFactory;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ReadStatMerger implements ApplicationCommand {

	
	private ApplicationOptions applicationOptions;
	
	private Logger logger;
	
	@Autowired
	public ReadStatMerger(LoggerService loggerService, ApplicationOptions applicationOptions, ReadPairStatisticsFactory readPairStatisticsFactory) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
	}
	
	@Override
	public void doWork() throws Exception {
		String fileOut = applicationOptions.getFileOut();
		File bamList = new File(applicationOptions.getBamList());
		
		if(!bamList.exists()) {
			throw new IllegalArgumentException("bamList "+ bamList.getAbsolutePath() + " not found!");
		}
		
		logger.info("Reading bam list: " + bamList.getAbsolutePath());
		List<String> lines = FileUtils.readLines(bamList);
		List<BAMInfo> infos = BAMInfo.getBAMInfoFromLines(lines);
	
		List<DelimitedFileEntry> entries = new LinkedList<DelimitedFileEntry>();
		
		logger.info("Parsing bam stats");
		for(BAMInfo info : infos) {
			String filePath = info.getBamStats();
			logger.info("Reading stat file for " + info.getSampleId());
			List<String> statLines = FileUtils.readLines(new File(filePath));
			
			if(statLines.size() != 2) {
				throw new RuntimeException("The stat file " +  filePath + " is not formatted correclty");
			}
			
			//line 1 should be header
			//line 2 should be the values
			
			DelimitedFileEntry dfe = DelimitedFileEntry.getDelimitedFileEntryFromDelimitedLine(info.getSampleId(),statLines.get(0),statLines.get(1),"\t");
			entries.add(dfe);
			
		}
		
		logger.info("Writing out merged matrix: " + fileOut);
		DelimitedFileEntry.writeToFile(entries, fileOut, "\t");
		
		
	}

}
