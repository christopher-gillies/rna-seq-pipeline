package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;

import org.springframework.stereotype.Component;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;


@Component
public class UniqueMappingFilter {

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
		
		while(iterator.hasNext()) {
			SAMRecord record = iterator.next();
			if(record.getMappingQuality() == 255) {
				writer.addAlignment(record);
			}
		}
		
	}
}
