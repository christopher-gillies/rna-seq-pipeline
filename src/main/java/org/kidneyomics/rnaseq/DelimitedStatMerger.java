package org.kidneyomics.rnaseq;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Scope("prototype")
@Component
public class DelimitedStatMerger implements ApplicationCommand {

	private BAMInfoService bamInfoService;
	
	private ApplicationOptions applicationOptions;
	
	private Logger logger;
	
	private String headerLineStartsWith = null;
	
	@Autowired
	public DelimitedStatMerger(LoggerService loggerService, ApplicationOptions applicationOptions, BAMInfoService bamInfoService) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
		this.bamInfoService = bamInfoService;
	}
	
	@Override
	public void doWork() throws Exception {
		String fileOut = applicationOptions.getFileOut();
		String bamList = applicationOptions.getBamList();
		

		logger.info("Reading bam list: " + bamList);
		List<BAMInfo> infos = bamInfoService.getBAMInfosFromBamList(bamList);
	
		List<DelimitedFileEntry> entries = new LinkedList<DelimitedFileEntry>();
		
		logger.info("Parsing bam stats");
		for(BAMInfo info : infos) {
			String filePath = info.getBamStats();
			logger.info("Reading stat file for " + info.getSampleId());
			List<String> statLines = FileUtils.readLines(new File(filePath));
			
			/*
			 * loop through the lines until we find a line that starts with headerLineStartsWith string
			 * break through the loop when it is found
			 */
			if(headerLineStartsWith != null) {
				Iterator<String> iter = statLines.iterator();
				while(iter.hasNext()) {
					String line = iter.next();
					if(line.startsWith(headerLineStartsWith)) {
						break;
					} else {
						iter.remove();
					}
				}
				
				if(statLines.size() == 0) {
					throw new RuntimeException("Could not find a header line starting with  " + headerLineStartsWith + " in " +  filePath );
				}
			}
			
			if(statLines.size() < 2) {
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

	String getHeaderLineStartsWith() {
		return headerLineStartsWith;
	}

	void setHeaderLineStartsWith(String headerLineStartsWith) {
		this.headerLineStartsWith = headerLineStartsWith;
	}
	
	

}
