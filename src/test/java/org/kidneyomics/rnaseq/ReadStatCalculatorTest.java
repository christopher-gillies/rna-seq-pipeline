package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.RandomUtils;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;
import org.kidneyomics.rnaseq.stats.ReadPairStatisticsFactory;
import org.slf4j.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMFileHeader.SortOrder;

public class ReadStatCalculatorTest {

	@Test
	public void test() throws Exception {
		LoggerService loggerService = new LoggerService();
		
		ApplicationOptions options = new ApplicationOptions(loggerService);
		
		ReadPairStatisticsFactory factory = new ReadPairStatisticsFactory();
		
		ReadStatCalculator rsc = new ReadStatCalculator(loggerService, options, factory);
		
		Logger logger = loggerService.getLogger(this);
		
		String tmpdir = FileUtils.getTempDirectoryPath();
		
		File sam = writeSam(logger, tmpdir);
		
		File log = new File(tmpdir + "/stat.log.txt");
		
		options.setFileIn(sam.getAbsolutePath());
		options.setFileOut(log.getAbsolutePath());
		
		rsc.doWork();
		
		assertTrue(log.exists());
		
		assertTrue(sam.exists());
		
		List<String> linesSam = FileUtils.readLines(sam);
		for(String line : linesSam) {
			System.err.println(line);
		}
		
		List<String> lines = FileUtils.readLines(log);
		assertEquals(2,lines.size());
		
		System.err.println(lines.get(0));
		System.err.println(lines.get(1));
		
		
		String[] header = lines.get(0).split("\t");
		
		assertEquals("A",header[0]);
		assertEquals("T",header[1]);
		assertEquals("G",header[2]);
		assertEquals("C",header[3]);
		assertEquals("N",header[4]);
		assertEquals("GC",header[5]);
		assertEquals("GT",header[6]);
		assertEquals("AT",header[7]);
		assertEquals("AC",header[8]);
		assertEquals("MEAN_PHRED_BASE_QUALITY",header[9]);
		assertEquals("MEAN_BASES_PER_READ_GREATER_THAN_Q30",header[10]);
		assertEquals("MEAN_INSERT_SIZE",header[11]);
		assertEquals("SD_INSERT_SIZE",header[12]);
		assertEquals("AAAAAAAAA",header[13]);
		
		
		String[] vals = lines.get(1).split("\t");
		
		assertEquals("0.5",vals[0]);
		assertEquals("0.5",vals[1]);
		assertEquals("0.0",vals[2]);
		assertEquals("0.0",vals[3]);
		assertEquals("0.0",vals[4]);
		assertEquals("0.0",vals[5]);
		assertEquals("0.5",vals[6]);
		assertEquals("1.0",vals[7]);
		assertEquals("0.5",vals[8]);
		assertEquals("39.0",vals[9]);
		assertEquals("78.0",vals[10]);
		assertEquals("158.15384615384616",vals[11]);
		assertEquals("19.61161351381841",vals[12]);
		assertEquals("1612.0",vals[13]);
		
		
		sam.delete();
		log.delete();
		
	}
	
	private File writeSam(Logger logger, String tmpdir) {
		File samOut = new File(tmpdir + "/test.sam");
		SAMFileHeader header = new SAMFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMSequenceRecord srec = new SAMSequenceRecord("1",10000);
		header.addSequence(srec);
		logger.info(srec.getSequenceIndex() + "");
		logger.info("Writing out to " + samOut.getAbsolutePath());
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileWriter writer = factory.makeSAMWriter(header, true, samOut);
		
		List<SAMRecord> records = new ArrayList<SAMRecord>(100);
		
		for(int i = 1; i <= 25; i++) {
			int pos = RandomUtils.nextInt(1, 100);
			{
				SAMRecord record = new SAMRecord(header);
				record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:" + i);
				record.setReadString(StringUtils.repeat("A", 39));
				record.setBaseQualityString(StringUtils.repeat("H", 39));
				record.setAlignmentStart(100 + pos);
				record.setMateAlignmentStart(100 + pos + 200);
				record.setCigarString("39M");
				record.setReferenceName("1");
				record.setReferenceIndex(0);
				record.setMateReferenceName("1");
				record.setMateReferenceIndex(0);
				record.setReadPairedFlag(true);
				record.setProperPairFlag(true);
				record.setReadUnmappedFlag(false);
				record.setMateUnmappedFlag(false);
				record.setReadNegativeStrandFlag(false);
				record.setMateNegativeStrandFlag(true);
				record.setFirstOfPairFlag(true);
				//logger.info(record.toString());
				//logger.info("Start: " + record.getStart());
				//logger.info("End: " + record.getAlignmentEnd());
				//writer.addAlignment(record);
				records.add(record);
			}
			{
				SAMRecord record = new SAMRecord(header);
				record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:" + i);
				record.setReadString(StringUtils.repeat("A", 39));
				record.setBaseQualityString(StringUtils.repeat("H", 39));
				record.setAlignmentStart(100 + pos + 200);
				record.setMateAlignmentStart(100 + pos);
				record.setCigarString("39M");
				record.setReferenceName("1");
				record.setReferenceIndex(0);
				record.setMateReferenceName("1");
				record.setMateReferenceIndex(0);
				record.setProperPairFlag(true);
				record.setReadPairedFlag(true);
				record.setReadUnmappedFlag(false);
				record.setMateUnmappedFlag(false);
				record.setReadNegativeStrandFlag(true);
				record.setMateNegativeStrandFlag(false);
				record.setSecondOfPairFlag(true);
				//logger.info(record.toString());
				//logger.info("Start: " + record.getStart());
				//logger.info("End: " + record.getAlignmentEnd());
				//writer.addAlignment(record);
				records.add(record);
			}
		}
		
		int pos = RandomUtils.nextInt(1, 100);
		{
			SAMRecord record = new SAMRecord(header);
			record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:" + 2222);
			record.setReadString(StringUtils.repeat("A", 39));
			record.setBaseQualityString(StringUtils.repeat("H", 39));
			record.setAlignmentStart(100 + pos);
			record.setMateAlignmentStart(100 + pos + 100);
			record.setCigarString("39M");
			record.setReferenceName("1");
			record.setReferenceIndex(0);
			record.setMateReferenceName("1");
			record.setMateReferenceIndex(0);
			record.setReadPairedFlag(true);
			record.setProperPairFlag(true);
			record.setReadUnmappedFlag(false);
			record.setMateUnmappedFlag(false);
			record.setReadNegativeStrandFlag(false);
			record.setMateNegativeStrandFlag(true);
			record.setFirstOfPairFlag(true);
			//logger.info(record.toString());
			//logger.info("Start: " + record.getStart());
			//logger.info("End: " + record.getAlignmentEnd());
			//writer.addAlignment(record);
			records.add(record);
		}
		{
			SAMRecord record = new SAMRecord(header);
			record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:" + 2222);
			record.setReadString(StringUtils.repeat("A", 39));
			record.setBaseQualityString(StringUtils.repeat("H", 39));
			record.setAlignmentStart(100 + pos + 100);
			record.setMateAlignmentStart(100 + pos);
			record.setCigarString("39M");
			record.setReferenceName("1");
			record.setReferenceIndex(0);
			record.setMateReferenceName("1");
			record.setMateReferenceIndex(0);
			record.setProperPairFlag(true);
			record.setReadPairedFlag(true);
			record.setReadUnmappedFlag(false);
			record.setMateUnmappedFlag(false);
			record.setReadNegativeStrandFlag(true);
			record.setMateNegativeStrandFlag(false);
			record.setSecondOfPairFlag(true);
			//logger.info(record.toString());
			//logger.info("Start: " + record.getStart());
			//logger.info("End: " + record.getAlignmentEnd());
			//writer.addAlignment(record);
			records.add(record);
		}

		Collections.sort(records, new SAMRecordCoordinateComparator());
		for(SAMRecord rec : records) {
			writer.addAlignment(rec);
		}
		

		writer.close();
		
		return samOut;
	}

}
