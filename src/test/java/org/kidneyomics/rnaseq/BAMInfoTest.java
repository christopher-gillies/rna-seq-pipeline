package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import org.junit.Test;

public class BAMInfoTest {

	@Test
	public void test() {
		String line = "25969	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//duplicate.output.metrics";
		BAMInfo b = BAMInfo.getBAMInfoFromLine(line);
		assertEquals("25969",b.sampleId);
		assertEquals("/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//final.bam",b.bam);
		assertEquals("/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//Log.final.out",b.log1);
		assertEquals("/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//Log.final.out",b.log2);
		assertEquals("/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//SJ.out.tab",b.splice1);
		assertEquals("/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//SJ.out.tab",b.splice2);
		assertEquals("/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//duplicate.output.metrics",b.dupmetrics);
	}

}
