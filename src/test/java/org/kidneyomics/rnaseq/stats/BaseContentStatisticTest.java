package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;
import org.kidneyomics.rnaseq.SAMRecordPair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BaseContentStatisticTest {

	@Test
	public void test() {
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r1 = new SAMRecord(header);
		r1.setReadString("AAATTTTGGGGGGCCCCC");
		
		SAMRecord r2 = new SAMRecord(header);
		r2.setReadString("AATTGGCCN");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.addPair(r1);
		pair.addPair(r2);
		
		
		
		BaseContentStatistic stat = new BaseContentStatistic();
		
		stat.addReadPair(pair);
		
		List<Double> result = stat.getStatistic();
		assertEquals(9,result.size());
		
		assertEquals("A\tT\tG\tC\tN\tGC\tGT\tAT\tAC",stat.getHeader());
		
		//A
		assertEquals(5 / (double) 27, result.get(0), 0.00001);
		//T
		assertEquals(6 / (double) 27, result.get(1), 0.00001);
		//G
		assertEquals(8 / (double) 27, result.get(2), 0.00001);
		//C
		assertEquals(7 / (double) 27, result.get(3), 0.00001);
		//N
		assertEquals(1 / (double) 27, result.get(4), 0.00001);
		
		//GC
		assertEquals(15 / (double) 27, result.get(5), 0.00001);
		//GT
		assertEquals(14 / (double) 27, result.get(6), 0.00001);
		//AT
		assertEquals(11 / (double) 27, result.get(7), 0.00001);
		//AC
		assertEquals(12 / (double) 27, result.get(8), 0.00001);
		
		
		//add second
		
		stat.addReadPair(pair);
		
		List<Double> result2 = stat.getStatistic();
		
		//A
		assertEquals(10 / (double) 54, result2.get(0), 0.00001);
		//T
		assertEquals(12 / (double) 54, result2.get(1), 0.00001);
		//G
		assertEquals(16 / (double) 54, result2.get(2), 0.00001);
		//C
		assertEquals(14 / (double) 54, result2.get(3), 0.00001);
		//N
		assertEquals(2 / (double) 54, result2.get(4), 0.00001);
		
		//GC
		assertEquals(30 / (double) 54, result2.get(5), 0.00001);
		//GT
		assertEquals(28 / (double) 54, result2.get(6), 0.00001);
		//AT
		assertEquals(22 / (double) 54, result2.get(7), 0.00001);
		//AC
		assertEquals(24 / (double) 54, result2.get(8), 0.00001);
		
		
		
	}

}
