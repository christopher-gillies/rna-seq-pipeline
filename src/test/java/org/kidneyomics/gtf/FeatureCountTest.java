package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;

public class FeatureCountTest {

	@Test
	public void test() {
		String line = "1	protein_coding	exon	183217372	183274009	.	-	.	transcript_id \"ENST00000294868\"; gene_id \"ENSG00000157064\";";
		Feature f = GTFFeatureBuilder.createFromLine(line, true);
		
		FeatureCount fc = new FeatureCount(f);
		
		assertEquals("ENSG00000157064_1_183217372_183274009", fc.getId());
		
		fc.addToCount(100);
		
		System.out.println("Length " + f.location().length());
		
		double length = 56638.0;
		double complexity = 2000000.0;
		double rpkm = 100.0 / length / complexity * Math.pow(10, 9);
		
		System.out.println("RPKM " + rpkm);
		
		assertEquals(rpkm,fc.getRPKM(complexity),0.0001);
	}

}
