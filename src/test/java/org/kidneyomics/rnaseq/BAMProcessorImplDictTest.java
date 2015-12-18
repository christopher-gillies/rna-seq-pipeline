package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.RandomUtils;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceRecord;

public class BAMProcessorImplDictTest {

	Logger logger = LoggerFactory.getLogger(BAMProcessorImplDictTest.class);
	
	@Test
	public void test() throws Exception {
		String tmpdir = FileUtils.getTempDirectoryPath();
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

		Collections.sort(records, new SAMRecordCoordinateComparator());
		for(SAMRecord rec : records) {
			writer.addAlignment(rec);
		}
		

		writer.close();
		
		
		String lines = FileUtils.readFileToString(samOut);
		
		System.out.println(lines);
		
		
		try(BAMProcessor processor = BAMProcessorImplDict.getBAMProcessor(samOut).withLogSkipSize(1)) {
			SAMRecordPair pair = null;
			int count = 0;
			while( ( pair = processor.getNextReadPair()) != null) {
			
				if(pair.getMate1() != null) {
					logger.info(pair.getMate1().getReadName());
					logger.info(pair.getMate1().getAlignmentStart() + "");
					count++;
				}
				if(pair.getMate2() != null) {
					logger.info(pair.getMate2().getReadName());
					logger.info(pair.getMate2().getAlignmentStart() + "");
					count++;
				}
			}
			
			assertEquals(50,count);
		}
		
		samOut.delete();
	}

	
	
	@Test
	public void test3() throws Exception {
		String tmpdir = FileUtils.getTempDirectoryPath();
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

		{
			SAMRecord record = new SAMRecord(header);
			record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:777");
			record.setReadString(StringUtils.repeat("A", 39));
			record.setBaseQualityString(StringUtils.repeat("H", 39));
			record.setAlignmentStart(100 + 77);
			record.setMateAlignmentStart(100 + 77 + 200);
			record.setCigarString("39M");
			record.setReferenceName("1");
			record.setReferenceIndex(0);
			record.setMateReferenceName("1");
			record.setMateReferenceIndex(0);
			record.setReadPairedFlag(true);
			record.setProperPairFlag(false);
			record.setReadUnmappedFlag(false);
			record.setMateUnmappedFlag(true);
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
			record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:777");
			record.setReadString(StringUtils.repeat("A", 39));
			record.setBaseQualityString(StringUtils.repeat("H", 39));
			record.setAlignmentStart(100 + 77 + 200);
			record.setMateAlignmentStart(100 + 77);
			record.setCigarString("39M");
			record.setReferenceName("1");
			record.setReferenceIndex(0);
			record.setMateReferenceName("1");
			record.setMateReferenceIndex(0);
			record.setProperPairFlag(false);
			record.setReadPairedFlag(true);
			record.setReadUnmappedFlag(true);
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
		
		
		String lines = FileUtils.readFileToString(samOut);
		
		System.out.println(lines);
		
		
		try(BAMProcessor processor = BAMProcessorImplDict.getBAMProcessor(samOut)) {
			SAMRecordPair pair = null;
			int count = 0;
			while( ( pair = processor.getNextReadPair()) != null) {
			
				if(pair.getMate1() != null) {
					logger.info(pair.getMate1().getReadName());
					logger.info(pair.getMate1().getAlignmentStart() + "");
					count++;
				}
				if(pair.getMate2() != null) {
					logger.info(pair.getMate2().getReadName());
					logger.info(pair.getMate2().getAlignmentStart() + "");
					count++;
				}
			}
			
			assertEquals(50,count);
		}
	}
	
	@Test
	public void test2() throws Exception {
		String tmpdir = FileUtils.getTempDirectoryPath();
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

		{
			SAMRecord record = new SAMRecord(header);
			record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:777");
			record.setReadString(StringUtils.repeat("A", 39));
			record.setBaseQualityString(StringUtils.repeat("H", 39));
			//record.setAlignmentStart(100 + 77);
			record.setMateAlignmentStart(100 + 77 + 200);
			//record.setCigarString("39M");
			//record.setReferenceName("1");
			//record.setReferenceIndex(0);
			record.setMateReferenceName("1");
			record.setMateReferenceIndex(0);
			
			
			record.setReadPairedFlag(true);
			record.setProperPairFlag(false);
			record.setReadUnmappedFlag(true);
			//record.setMateUnmappedFlag(false);
			//record.setReadNegativeStrandFlag(false);
			//record.setMateNegativeStrandFlag(true);
			record.setFirstOfPairFlag(true);
			//logger.info(record.toString());
			//logger.info("Start: " + record.getStart());
			//logger.info("End: " + record.getAlignmentEnd());
			//writer.addAlignment(record);
			records.add(record);
		}
		{
			SAMRecord record = new SAMRecord(header);
			record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:777");
			record.setReadString(StringUtils.repeat("A", 39));
			record.setBaseQualityString(StringUtils.repeat("H", 39));
			record.setAlignmentStart(100 + 77 + 200);
			//record.setMateAlignmentStart(100 + 77);
			record.setCigarString("39M");
			record.setReferenceName("1");
			record.setReferenceIndex(0);
			//record.setMateReferenceName("1");
			//record.setMateReferenceIndex(0);
			record.setProperPairFlag(false);
			record.setReadPairedFlag(true);
			//record.setReadUnmappedFlag(true);
			record.setMateUnmappedFlag(true);
			record.setReadNegativeStrandFlag(true);
			//record.setMateNegativeStrandFlag(false);
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
		
		
		String lines = FileUtils.readFileToString(samOut);
		
		System.out.println(lines);
		
		
		try(BAMProcessor processor = BAMProcessorImplDict.getBAMProcessor(samOut)) {
			SAMRecordPair pair = null;
			int count = 0;
			while( ( pair = processor.getNextReadPair()) != null) {
			
				if(pair.getMate1() != null) {
					logger.info(pair.getMate1().getReadName());
					logger.info(pair.getMate1().getAlignmentStart() + "");
					count++;
				}
				if(pair.getMate2() != null) {
					logger.info(pair.getMate2().getReadName());
					logger.info(pair.getMate2().getAlignmentStart() + "");
					count++;
				}
			}
			
			assertEquals(50,count);
		}
	}
}
