package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import static org.mockito.Mockito.*;

public class TestUniqueMappingFilter {

	@Test
	public void test() {
		
		LoggerService loggerService = mock(LoggerService.class);
		Logger logger = LoggerFactory.getLogger(TestUniqueMappingFilter.class);
		when(loggerService.getLogger(anyObject())).thenReturn(logger);
		
		UniqueMappingFilter filter = new UniqueMappingFilter(loggerService);
		
		SAMRecordIterator iterator = mock(SAMRecordIterator.class);
		
		SAMFileWriter writer = mock(SAMFileWriter.class);
		
		
				
		when(iterator.hasNext()).thenReturn(true, true, true, false);
		
		SAMRecord r1 = mock(SAMRecord.class);
		when(r1.getMappingQuality()).thenReturn(255);
		
		
		SAMRecord r2 = mock(SAMRecord.class);
		when(r2.getMappingQuality()).thenReturn(1);
		
		SAMRecord r3 = mock(SAMRecord.class);
		when(r3.getMappingQuality()).thenReturn(255);
		
		when(iterator.next()).thenReturn(r1,r2,r3);
		
		filter.filterIterator(iterator, writer);

		verify(writer).addAlignment(r1);
		verify(writer).addAlignment(r3);
		verify(writer,atMost(2)).addAlignment((SAMRecord) anyObject());
		
	}

}
