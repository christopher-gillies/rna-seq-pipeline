package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import org.junit.Test;
import org.kidneyomics.rnaseq.SAMRecordPair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class MeanBasesPerReadGreaterThanQ30Test {

	@Test
	public void test() {
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r1 = new SAMRecord(header);
		r1.setBaseQualityString("??@;??DAFHHH");
		
		SAMRecord r2 = new SAMRecord(header);
		r2.setBaseQualityString("?;=;DDAD?");
		
		/*
		 * = = 61
		 * ; = 59
		 * ? = 63
		 * @ = 64
		 * D = 68
		 * F = 70
		 * H = 72
		 * A = 65
		 */
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.addPair(r1);
		pair.addPair(r2);
		
		MeanBasesPerReadGreaterThanQ30 stat = new MeanBasesPerReadGreaterThanQ30();
		
		stat.addReadPair(pair);
		
		
		assertEquals(1,stat.getStatistic().size());
		
		System.err.println(stat.getStatistic().get(0));
		assertEquals(11.0, stat.getStatistic().get(0), 0.0001);
		
		assertEquals("MEAN_BASES_PER_READ_GREATER_THAN_Q30",stat.getHeader());
	}

	
	@Test
	public void test2() {
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r1 = new SAMRecord(header);
		r1.setBaseQualityString("??@;??DAFHHH");
		
		SAMRecord r2 = new SAMRecord(header);
		r2.setBaseQualityString("?;=;DDAD?");
		
		/*
		 * = = 61
		 * ; = 59
		 * ? = 63
		 * @ = 64
		 * D = 68
		 * F = 70
		 * H = 72
		 * A = 65
		 */
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.addPair(r1);
		pair.addPair(r2);
		
		MeanBasesPerReadGreaterThanQ30 stat = new MeanBasesPerReadGreaterThanQ30();
		
		stat.addReadPair(pair);
		
		SAMRecord r3 = new SAMRecord(header);
		r3.setBaseQualityString("??@;??DAFHHH");
		
		SAMRecord r4 = new SAMRecord(header);
		r4.setBaseQualityString("?;=;?DAD?");
		
		SAMRecordPair pair2 = new SAMRecordPair();
		pair2.addPair(r3);
		pair2.addPair(r4);
		
		stat.addReadPair(pair2);
		
		assertEquals(1,stat.getStatistic().size());
		
		System.err.println(stat.getStatistic().get(0));
		assertEquals(10.5, stat.getStatistic().get(0), 0.0001);
	}
	
}
