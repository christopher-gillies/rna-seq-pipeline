package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;
import static org.mockito.Mockito.*;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class DelimitedStatMergerTest {

	@Test
	public void test() throws Exception {
		LoggerService ls = new LoggerService();
		
		Logger logger = ls.getLogger(this);
		ApplicationOptions opts = new ApplicationOptions(ls);
		
		
		File tmpDir = FileUtils.getTempDirectory();
		
		File out = new File(tmpDir.getAbsolutePath() + "/test.stats");
		
		
		
		opts.setFileOut(out.getAbsolutePath());
		Resource r = new ClassPathResource("duplicate.metrics.test");
		
		String res = FileUtils.readFileToString(r.getFile());
		
		System.err.println(res);
		
		BAMInfoService bs = mock(BAMInfoService.class);
		
		List<BAMInfo> list = new LinkedList<BAMInfo>();
		
		BAMInfo info1 = new BAMInfo();
		
		info1.setSampleId("sample1");
		//info1.setBamStats(r.getFile().getAbsolutePath());
		
		info1.setDupmetrics(r.getFile().getAbsolutePath());
		
		logger.info(info1.getBamStats());
		
		list.add(info1);
		list.add(info1);
		
		when(bs.getBAMInfosFromBamList(anyString())).thenReturn(list);
		
		DelimitedStatMerger merger = new DelimitedStatMerger(ls, opts, bs);
		merger.setHeaderLineStartsWith("LIBRARY");
		merger.setBamInfoSelector(new BAMInfoDuplicateSelector());
		
		merger.doWork();
		
		assertTrue(out.exists());
		
		List<String> lines = FileUtils.readLines(out);
		
		assertEquals(3,lines.size());
		
		assertEquals("SAMPLE_ID	LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE",lines.get(0));
		assertEquals("sample1	26003	0	25838752	0	0	14531382	7050920	0.56238	16793725",lines.get(1));
		assertEquals("sample1	26003	0	25838752	0	0	14531382	7050920	0.56238	16793725",lines.get(2));
		
		out.delete();
	}

}
