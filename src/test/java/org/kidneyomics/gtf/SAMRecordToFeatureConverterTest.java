package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMFileHeader.SortOrder;

public class SAMRecordToFeatureConverterTest {

	static SAMRecordToFeatureConverter converter = new SAMRecordToFeatureConverter();
	@Test
	public void test() {
		
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMSequenceRecord srec = new SAMSequenceRecord("1",10000);
		header.addSequence(srec);
		
		SAMRecord record = new SAMRecord(header);
		record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:1");
		record.setReadString(StringUtils.repeat("A", 39));
		record.setBaseQualityString(StringUtils.repeat("H", 39));
		record.setAlignmentStart(100);
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
		
		
		List<Feature> features = converter.convert(record);
		
		assertTrue(features.size() == 1);
		
		assertEquals(100,features.get(0).location().bioStart());
		
		assertEquals(138,features.get(0).location().bioEnd());
		
		assertEquals('+',features.get(0).location().bioStrand());
		
		assertEquals(record.getAlignmentEnd(),features.get(0).location().bioEnd());
		
	}
	
	
	@Test
	public void test2() {
		
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMSequenceRecord srec = new SAMSequenceRecord("1",10000);
		header.addSequence(srec);
		
		SAMRecord record = new SAMRecord(header);
		record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:1");
		record.setReadString(StringUtils.repeat("A", 39));
		record.setBaseQualityString(StringUtils.repeat("H", 39));
		record.setAlignmentStart(100);
		record.setCigarString("3S36M");
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
		
		

		List<Feature> features = converter.convert(record);
		
		assertTrue(features.size() == 1);
		
		assertEquals(100,features.get(0).location().bioStart());
		
		assertEquals(135,features.get(0).location().bioEnd());
		
		assertEquals('+',features.get(0).location().bioStrand());
		
		assertEquals(record.getAlignmentEnd(),features.get(0).location().bioEnd());
		
	}
	
	
	@Test
	public void test3() {
		
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMSequenceRecord srec = new SAMSequenceRecord("1",10000);
		header.addSequence(srec);
		
		SAMRecord record = new SAMRecord(header);
		record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:1");
		record.setReadString(StringUtils.repeat("A", 39));
		record.setBaseQualityString(StringUtils.repeat("H", 39));
		record.setAlignmentStart(100);
		record.setCigarString("10M5N10M");
		record.setReferenceName("1");
		record.setReferenceIndex(0);
		record.setMateReferenceName("1");
		record.setMateReferenceIndex(0);
		record.setReadPairedFlag(true);
		record.setProperPairFlag(true);
		record.setReadUnmappedFlag(false);
		record.setMateUnmappedFlag(false);
		record.setReadNegativeStrandFlag(true);
		record.setMateNegativeStrandFlag(false);
		record.setFirstOfPairFlag(true);
		
		

		List<Feature> features = converter.convert(record);
		
		assertTrue(features.size() == 2);
		
		assertEquals(100,features.get(0).location().bioStart());
		
		assertEquals(109,features.get(0).location().bioEnd());
		
		assertEquals('-',features.get(0).location().bioStrand());
		
		
		assertEquals(115,features.get(1).location().bioStart());
		
		assertEquals(124,features.get(1).location().bioEnd());
		
		assertEquals('-',features.get(1).location().bioStrand());
		
		
		assertEquals(record.getAlignmentEnd(),features.get(1).location().bioEnd());
		
	}
	
	
	@Test
	public void test4() {
		
		
		SAMFileHeader header = new SAMFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMSequenceRecord srec = new SAMSequenceRecord("1",10000);
		header.addSequence(srec);
		
		SAMRecord record = new SAMRecord(header);
		record.setReadName("D7DHSVN1:225:D2HR7ACXX:4:1107:16208:1");
		record.setReadString(StringUtils.repeat("A", 39));
		record.setBaseQualityString(StringUtils.repeat("H", 39));
		record.setAlignmentStart(100);
		record.setCigarString("10M500N10M500N10M");
		record.setReferenceName("1");
		record.setReferenceIndex(0);
		record.setMateReferenceName("1");
		record.setMateReferenceIndex(0);
		record.setReadPairedFlag(true);
		record.setProperPairFlag(true);
		record.setReadUnmappedFlag(false);
		record.setMateUnmappedFlag(false);
		record.setReadNegativeStrandFlag(true);
		record.setMateNegativeStrandFlag(false);
		record.setFirstOfPairFlag(true);
		
		

		List<Feature> features = converter.convert(record);
		
		assertTrue(features.size() == 3);
		
		assertEquals(100,features.get(0).location().bioStart());
		
		assertEquals(109,features.get(0).location().bioEnd());
		
		assertEquals('-',features.get(0).location().bioStrand());
		
		
		assertEquals(610,features.get(1).location().bioStart());
		
		assertEquals(619,features.get(1).location().bioEnd());
		
		assertEquals('-',features.get(1).location().bioStrand());
		
		
		assertEquals(1120,features.get(2).location().bioStart());
		
		assertEquals(1129,features.get(2).location().bioEnd());
		
		assertEquals('-',features.get(2).location().bioStrand());
		
		
		assertEquals(record.getAlignmentEnd(),features.get(2).location().bioEnd());
		
	}


}
