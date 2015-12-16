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
	
	
	@Test
	public void test3() {
		
		logger.info("BiojavaFeatureTest#test3");
		Location l1 = Location.fromBio(100, 200, '+');
		Location l2 = Location.fromBio(1, 99, '+');
		
		Location intersection = l1.intersection(l2);
		
		logger.info(intersection.length() + "");
		//assertNull(intersection);
		
		Location union = l1.union(l2);
		
		logger.info(union.bioStart() + "");
		logger.info(union.bioEnd() + "");
		logger.info(union.length() + "");
	}
	
	@Test
	public void test4() {
		
		logger.info("BiojavaFeatureTest#test4");
		Location l1 = Location.fromBio(100, 200, '+');
		Location l2 = Location.fromBio(1, 99, '+');
		
		assertFalse(l1.overlaps(l2));
		assertTrue(l1.isAfter(l2));
		assertTrue(l2.isBefore(l1));
	}
	
	@Test
	public void test5() {
		
		logger.info("BiojavaFeatureTest#test5");
		Location l1 = Location.fromBio(100, 200, '+');
		Location l2 = Location.fromBio(1, 100, '+');
		
		assertTrue(l1.overlaps(l2));
		assertFalse(l2.isBefore(l1));
		assertFalse(l1.isAfter(l2));
	}
	
	@Test
	public void test6() {
		
		logger.info("BiojavaFeatureTest#test6");
		Location l1 = Location.fromBio(51227320, 51227381, '+');
		Location l2 = Location.fromBio(51227323, 51227382, '+');
		
		Location union = l1.union(l2);
		assertEquals(51227320,union.bioStart());
		assertEquals(51227382,union.bioEnd());
	}
	
	
	@Test
	public void test6_5() {
		
		logger.info("BiojavaFeatureTest#test6_5");
		Location l1 = Location.fromBio(99, 100, '+');
		Location l2 = Location.fromBio(99, 101, '+');
		
		Location union = l1.union(l2);
		assertEquals(99,union.bioStart());
		assertEquals(101,union.bioEnd());
	}
	
	//@Test
	public void test7() {
		
		logger.info("BiojavaFeatureTest#test7");
		Location l1 = Location.fromBio(100, 200, '+');
		Location l2 = Location.fromBio(1, 99, '+');
		Location intersection = l1.intersection(l2);
		assertNull(intersection);
	}
	
	
	@Test
	public void test8() {
		
		logger.info("BiojavaFeatureTest#test6");
		Location l1 = Location.fromBio(51227320, 51227380, '+');
		Location l2 = Location.fromBio(51227323, 51227382, '+');
		
		Location union = l1.union(l2);
		assertEquals(51227320,union.bioStart());
		assertEquals(51227382,union.bioEnd());
	}
	

}
