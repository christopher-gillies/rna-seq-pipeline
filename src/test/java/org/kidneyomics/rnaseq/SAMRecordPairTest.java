package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class SAMRecordPairTest {

	@Test
	public void test() {
		SAMFileHeader header = new SAMFileHeader();
		
		SAMRecordPair pair = new SAMRecordPair();
		
		SAMRecord r1 = new SAMRecord(header);
		r1.setAlignmentStart(100);
		
		SAMRecord r2 = new SAMRecord(header);
		r2.setAlignmentStart(150);
		
		pair.addPair(r1);
		pair.addPair(r2);
		
		pair.reorderMatesByCoordinate();
		
		assertTrue(pair.getMate1() == r1);
		
		assertTrue(pair.getMate2() == r2);
		
		r1.setAlignmentStart(300);
		
		pair.reorderMatesByCoordinate();
		
		assertTrue(pair.getMate1() == r2);
		
		assertTrue(pair.getMate2() == r1);
		
		r1.setAlignmentStart(150);
		
		pair.reorderMatesByCoordinate();
		
		assertTrue(pair.getMate1() == r2);
		
		assertTrue(pair.getMate2() == r1);
		
	}

}
