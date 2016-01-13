package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import org.junit.Test;
import org.kidneyomics.rnaseq.SAMRecordPair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class InsertSizeStatisticTest {

	@Test
	public void test() {
		
		
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r1 = new SAMRecord(header);
		
		r1.setAlignmentStart(100);
		r1.setCigarString("51M");
		

		
		SAMRecord r2 = new SAMRecord(header);
		
		r2.setAlignmentStart(30);
		r2.setCigarString("51M");
		
		
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.addPair(r1);
		pair.addPair(r2);
		
		SAMRecordPair pair2 = new SAMRecordPair();
		
		SAMRecord r3 = new SAMRecord(header);
		
		r3.setAlignmentStart(100);
		r3.setCigarString("51M");
		
		SAMRecord r4 = new SAMRecord(header);
		
		r4.setAlignmentStart(40);
		r4.setCigarString("51M");
		
		pair2.addPair(r3);
		pair2.addPair(r4);
		
		
		InsertSizeStatistic stat = new InsertSizeStatistic();
		
		stat.addReadPair(pair);
		stat.addReadPair(pair2);
		
		
		assertEquals(15, stat.getStatisticAsMap().get("MEAN_INSERT_SIZE"),0.001);
		assertEquals( 7.071068, stat.getStatisticAsMap().get("SD_INSERT_SIZE"),0.001);
	}

}
