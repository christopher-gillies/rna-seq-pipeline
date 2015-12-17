package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class MaxAlignmentDistanceFilterTest {

	@Test
	public void test1() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("nM", 5);
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setAttribute("nM", 1);
		
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		assertTrue(filter.keep(pair));
	}
	
	
	@Test
	public void test2() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("nM", 5);
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setAttribute("nM", 7);
		
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		assertTrue(filter.keep(pair));
	}
	
	@Test
	public void test3() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("nM", 5);
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setAttribute("nM", 8);
		
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		assertFalse(filter.keep(pair));
	}
	
	@Test
	public void test4() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("nM", 8);
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setAttribute("nM", 8);
		
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		assertFalse(filter.keep(pair));
	}

	@Test
	public void test5() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("nM", 8);
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setAttribute("nM", 0);
		
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		assertFalse(filter.keep(pair));
	}
	
	@Test
	public void test6() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("nM", 8);
		
		pair.setMate1(record1);
		
		assertFalse(filter.keep(pair));
	}
	
	@Test
	public void test7() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		record1.setAttribute("NM", 8);
		
		pair.setMate1(record1);
		
		assertFalse(filter.keep(pair));
	}
	
	@Test
	public void test8() {
		MaxAlignmentDistanceFilter filter = new MaxAlignmentDistanceFilter(7);
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		
		pair.setMate1(record1);
		
		assertTrue(filter.keep(pair));
	}
}
