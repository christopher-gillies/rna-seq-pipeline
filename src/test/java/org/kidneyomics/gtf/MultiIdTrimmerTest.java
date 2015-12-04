package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.junit.Test;

public class MultiIdTrimmerTest {

	@Test
	public void test() {
		String id = "ENST00000416570_ENST00000448975";
		String result = MultiIdTrimmer.trim(id);
		
		assertEquals("ENST00000416570",result);
	}
	
	@Test
	public void test2() {
		String ex = "1	lincRNA	transcript	763079	788146	.	+	.	transcript_id \"ENST00000416570_ENST00000448975\"; locus_id \"1:762988-794826W\"; gene_id \"ENSG00000228794\"; reads 0.000000; length 612; RPKM 0.000000";
		String result = MultiIdTrimmer.trim(ex);
		assertEquals("1	lincRNA	transcript	763079	788146	.	+	.	transcript_id \"ENST00000416570\"; locus_id \"1:762988-794826W\"; gene_id \"ENSG00000228794\"; reads 0.000000; length 612; RPKM 0.000000",result);
	}

}
