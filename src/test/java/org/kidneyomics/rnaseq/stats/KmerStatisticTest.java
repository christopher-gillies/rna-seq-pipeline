package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import java.util.Map;

import org.junit.Test;
import org.kidneyomics.rnaseq.SAMRecordPair;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class KmerStatisticTest {

	@Test
	public void test() {
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r1 = new SAMRecord(header);
		r1.setReadString("AAATTTTGGGGGGCCCCC");
		
		//AA AT TG GC CN TT GG CC
		SAMRecord r2 = new SAMRecord(header);
		r2.setReadString("AATTGGCCN");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.addPair(r1);
		pair.addPair(r2);
		
		KmerStatistic stat = new KmerStatistic(2);
		
		stat.addReadPair(pair);
		
		Map<String,Double> map = stat.getStatisticAsMap();
		
		assertEquals(8, map.size());
		
		assertEquals(3.0,map.get("AA"),0.0000000000001);
		assertEquals(2.0,map.get("AT"),0.0000000000001);
		assertEquals(2.0,map.get("TG"),0.0000000000001);
		assertEquals(2.0,map.get("GC"),0.0000000000001);
		assertEquals(1.0,map.get("CN"),0.0000000000001);
		
		assertEquals(4.0,map.get("TT"),0.0000000000001);
		assertEquals(6.0,map.get("GG"),0.0000000000001);
		assertEquals(5.0,map.get("CC"),0.0000000000001);
	}

}
