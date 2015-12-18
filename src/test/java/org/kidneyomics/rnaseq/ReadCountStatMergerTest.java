package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.context.support.ReloadableResourceBundleMessageSource;

public class ReadCountStatMergerTest {

	
	Logger logger = LoggerFactory.getLogger(ReadCountStatMergerTest.class);
	
	
	@Test
	public void test() throws IOException {
		String tmpDir = FileUtils.getTempDirectoryPath();
		
		Resource r1 = new ClassPathResource("test.1.gtf.stats");
		File gtf1Ref = r1.getFile();
		
		File gtf1Out = new File(tmpDir + "/test.1.gtf.stats");
		
		Resource r2 = new ClassPathResource("test.2.gtf.stats");
		File gtf2Ref = r2.getFile();
		File gtf2Out = new File(tmpDir + "/test.2.gtf.stats");
		
		
		FileUtils.copyFile(gtf1Ref, gtf1Out);
		FileUtils.copyFile(gtf2Ref, gtf2Out);
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("1\t" + gtf1Out.getAbsolutePath() + "\n");
		sb.append("2\t" + gtf2Out.getAbsolutePath() + "\n");
		
		File statList = new File(tmpDir + "/gtf.list" );
		
		File outfile = new File(tmpDir + "/merged.txt");
		
		FileUtils.write(statList, sb.toString());
		
		ApplicationOptions ao = new ApplicationOptions(new LoggerService());
		
		ao.setFileIn(statList.getAbsolutePath());
		ao.setFileOut(outfile.getAbsolutePath());
		
		logger.info("Out file: " + outfile.getAbsolutePath());
		ReadCountStatMerger readCountStatMerger = new ReadCountStatMerger(new LoggerService(), ao);
		readCountStatMerger.mergeStatFiles();
		
		
		assertTrue(statList.exists());
		List<String> lines = FileUtils.readLines(outfile);
		
		assertEquals(3,lines.size());
		
		String line1 = "#ID	TOTAL_READS	AMBIGUOUS_READS	AMBIGUOUS_FRACTION	FILTERED_OUT_READS	FILTERED_OUT_FRACTION	UNMAPPED_READS	UNMAPPED_FRACTION	PARTIALLY_MAPPED_READS	PARTIALLY_MAPPED_FRACTION	MAPPED_READS	MAPPED_FRACTION";
		assertEquals(line1,lines.get(0));
		
		String line2 = "1	19424775	140795	0.007248217804324632	54641	0.002812954075401131	4675060.45249581	0.2406751405097773	1135951	0.05847949332746454	14554278.547491638	0.7492636876098507";
		assertEquals(line2,lines.get(1));
		
		String line3 = "2	21459362	149389	0.006961483757066031	72382	0.003372980054113445	5090658.268649271	0.23722318811944507	1011788	0.047149025213331135	16146932.731341528	0.7524423480689467";
		assertEquals(line3,lines.get(2));
	}

}
