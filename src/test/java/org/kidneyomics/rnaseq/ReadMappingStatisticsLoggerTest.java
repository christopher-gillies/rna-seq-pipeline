package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import static org.mockito.Mockito.*;

public class ReadMappingStatisticsLoggerTest {

	@Test
	public void test() throws IOException {
		File f = new File( FileUtils.getTempDirectoryPath() + "/stats.log");
		
		FeatureCounter fc = mock(FeatureCounter.class);
		
		/*
		 * 
		 * 		double ambiguousReadCount = fc.getAmbiguousReadCount();
		double mappedReadCount = fc.getMappedReadCount();
		double totalReads = fc.getTotalCount();
		long partiallyMappedReads = fc.getNumberOfPartiallyUnmappedReads();
		double unmappedReadCount = fc.getUnmappedReadCount();
		 */
		
		when(fc.getAmbiguousReadCount()).thenReturn(37737.0);
		when(fc.getMappedReadCount()).thenReturn(1.6170759914677566E7);
		when(fc.getTotalCount()).thenReturn(2.1117682E7);
		when(fc.getNumberOfPartiallyUnmappedReads()).thenReturn(1873854l);
		when(fc.getUnmappedReadCount()).thenReturn(4909185.085325402);
		ReadMappingStatisticsLogger.writeStats(f, fc);
		
		List<String> lines = FileUtils.readLines(f);
		
		assertEquals(2,lines.size());
		assertFalse(lines.get(1).contains("E"));
		
	}

}
