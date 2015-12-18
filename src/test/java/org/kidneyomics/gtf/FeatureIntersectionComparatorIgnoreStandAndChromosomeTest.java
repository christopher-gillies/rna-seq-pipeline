package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;

public class FeatureIntersectionComparatorIgnoreStandAndChromosomeTest {

	@Test
	public void test() {
		FeatureIntersectionComparatorIgnoreStrandAndChromosome comparator = new FeatureIntersectionComparatorIgnoreStrandAndChromosome();
		
		Location l1 = Location.fromBio(100, 500, '+');
		Location l2 = Location.fromBio(300, 400, '-');
		Location l3 = Location.fromBio(200, 299, '+');
		
		Feature o1 = new Feature("chr1", "", "", l1, -1.0, -1, "");
		Feature o2 = new Feature("chr1", "", "", l2, -1.0, -1, "");
		Feature o3 = new Feature("chr1", "", "", l3, -1.0, -1, "");
		
		
		
		assertEquals(0,comparator.compare(o1, o2));
		assertEquals(1,comparator.compare(o2, o3));
		assertEquals(0,comparator.compare(o1, o3));
		
		
		assertEquals(-1,comparator.compare(o3, o2));
		assertEquals(0,comparator.compare(o2, o1));
		assertEquals(0,comparator.compare(o3, o1));
		
		assertEquals(0,comparator.compare(o1, o1));
		assertEquals(0,comparator.compare(o2, o2));
		assertEquals(0,comparator.compare(o3, o3));
		
		
		
	}

}
