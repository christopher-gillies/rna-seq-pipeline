package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;

import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;


@Component
public class UniqueMappingFilter implements ApplicationCommand {

	ApplicationOptions applicationOptions;
	
	Logger logger;
	
	@Autowired
	public UniqueMappingFilter(LoggerService loggerService, ApplicationOptions applicationOptions) {
		logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
	}
	
	public void filter(File in, File out) throws IOException {
		SamReader reader = SamReaderFactory.makeDefault().open(in);
		
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		
		SAMFileWriter writer = factory.makeBAMWriter(reader.getFileHeader(), false, out);
		
		SAMRecordIterator iterator = reader.iterator();

		filterIterator(iterator,writer);
		
		writer.close();
		
		reader.close();
		
	}
	
	
	public void filterIterator(SAMRecordIterator iterator, SAMFileWriter writer) {
		
		int counter = 0;
		int uniqueReads = 0;
		while(iterator.hasNext()) {
			counter++;
			if(counter % 1000000 == 0) {
				logger.info("Scanned " + counter + " records");
			} 
			SAMRecord record = iterator.next();
			if(record.getMappingQuality() == 255) {
				uniqueReads++;
				if(uniqueReads % 1000000 == 0) {
					logger.info(uniqueReads + " uniquely mapped records written of " + counter + " scanned");
				} 
				writer.addAlignment(record);
			}
		}
		
		logger.info("Found " + uniqueReads + " uniquely mapped records");
		
	}

	@Override
	public void doWork() throws Exception {
    	File in = new File(applicationOptions.getFileIn());
    	File out = new File(applicationOptions.getFileOut());
    	this.filter(in,out);
	}
}
