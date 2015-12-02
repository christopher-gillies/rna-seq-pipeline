package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BiojavaFeatureTest {

	Logger logger = LoggerFactory.getLogger(BiojavaFeatureTest.class);
	@Test
	public void test() {
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '-'), 0.0, 0, "");
		
		Location intersection = f1.location().intersection(f2.location().opposite());
		
		assertEquals(intersection.bioStart(),100);
		assertEquals(intersection.bioEnd(),150);
		assertEquals(intersection.length(), 150 - 100 + 1);
		
	}
	
	
	@Test
	public void test2() {
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		Location intersection = f1.location().intersection(f2.location());
		
		logger.info(intersection.bioStart() + "");
		logger.info(intersection.bioEnd() + "");
		logger.info(intersection.bioStrand() + "");
		
		assertEquals(intersection.bioStart(),100);
		assertEquals(intersection.bioEnd(),150);
		assertEquals(intersection.length(), 150 - 100 + 1);
	}

}
