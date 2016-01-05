package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import org.junit.Test;
import org.kidneyomics.rnaseq.SAMRecordPair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class MeanBaseQualityStatisticTest {

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
		
		MeanBaseQualityStatistic stat = new MeanBaseQualityStatistic();
		
		stat.addReadPair(pair);
		
		int count = 21;
		// "??@;??DAFHHH"
		int qualities = 63 + 63  + 64 + 59 + 63 + 63 + 68 + 65 + 70 + 72 + 72 + 72;
		// ?;=;DDAD?
		qualities += 63 + 59 + 61 + 59 + 68 + 68 + 65 + 68 + 63; 
		
		assertEquals(1,stat.getStatistic().size());
		
		System.err.println(stat.getStatistic().get(0));
		assertEquals(qualities / (double) count - 33, stat.getStatistic().get(0), 0.0001);
		
		assertEquals("MEAN_PHRED_BASE_QUALITY",stat.getHeader());
	}

}
