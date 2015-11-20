package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import org.junit.Test;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import static org.mockito.Mockito.*;

public class TestUniqueMappingFilter {

	@Test
	public void test() {
		UniqueMappingFilter filter = new UniqueMappingFilter();
		
		SAMRecordIterator iterator = mock(SAMRecordIterator.class);
		
		SAMFileWriter writer = mock(SAMFileWriter.class);
		
				
		when(iterator.hasNext()).thenReturn(true, true, false);
		
		SAMRecord r1 = mock(SAMRecord.class);
		when(r1.getMappingQuality()).thenReturn(255);
		
		
		SAMRecord r2 = mock(SAMRecord.class);
		when(r2.getMappingQuality()).thenReturn(1);
		
		
		when(iterator.next()).thenReturn(r1,r2);
		
		filter.filterIterator(iterator, writer);

		verify(writer).addAlignment(r1);
		verify(writer,atMost(1)).addAlignment((SAMRecord) anyObject());
		
	}

}
