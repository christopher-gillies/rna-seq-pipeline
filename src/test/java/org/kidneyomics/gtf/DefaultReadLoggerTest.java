package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;

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
		
		readLogger.logRead("test_1", record);
		
		readLogger.close();
		
		List<String> lines = FileUtils.readLines(f);
		for(String line : lines) {
			logger.info(line);
		}
		
		assertTrue(lines.size() == 2);
	
		f.delete();
	}

}
