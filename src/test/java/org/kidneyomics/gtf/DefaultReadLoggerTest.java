package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;

public class DefaultReadLoggerTest {

	Logger logger = LoggerFactory.getLogger(DefaultReadLoggerTest.class);
	
	@Test
	public void test() throws IOException {
		/*
		 * 			writer.write(type);
			writer.write("\t");
			writer.write(record.getReadName());
			writer.write("\t");
			writer.write(record.getReferenceName());
			writer.write("\t");
			writer.write(record.getAlignmentStart());
			writer.write("\t");
			writer.write(record.getCigarString());
			writer.write("\t");
			writer.write(record.getAlignmentEnd());
			writer.write("\t");
			writer.write(record.getReadLength());
			writer.write("\n");
		 */
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		record.setReadName("NAME");
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		logger.info("" + record.getAlignmentStart());
		record.setCigarString("39M");
		record.setReadString(StringUtil.repeatCharNTimes('A', 39));
		record.setBaseQualityString(StringUtil.repeatCharNTimes('H', 39));
		
		DefaultReadLogger readLogger = new DefaultReadLogger();
		
		File f = new File(FileUtils.getTempDirectoryPath() + "/log.test");
		
		readLogger.setFile(f);
		
		Map<Feature,Integer> mappedFeatures = new HashMap<>();
		
		Feature feature = new Feature("chr","","",Location.fromBio(100, 100, '+'),0.0,0,"id \"EXON_ID\"");
		mappedFeatures.put(feature, 13);
		
		readLogger.logRead("test_1", record,mappedFeatures);
		
		readLogger.close();
		
		List<String> lines = FileUtils.readLines(f);
		for(String line : lines) {
			logger.info(line);
		}
		
		assertTrue(lines.size() == 2);
	
		f.delete();
	}

}
