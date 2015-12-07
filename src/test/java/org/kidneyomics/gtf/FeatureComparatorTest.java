package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Location;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;

public class FeatureComparatorTest {

	@Test
	public void test() {
		FeatureComparator comparator = new FeatureComparator();
		
		Location l1 = Location.fromBio(100, 200, '+');
		Location l2 = Location.fromBio(150, 200, '-');
		Location l3 = Location.fromBio(160, 200, '-');
		
		Feature o1 = new Feature("chr1", "", "", l1, -1.0, -1, "");
		Feature o2 = new Feature("chr1", "", "", l2, -1.0, -1, "");
		Feature o3 = new Feature("chr1", "", "", l3, -1.0, -1, "");
		
		
		
		assertEquals(-1,comparator.compare(o1, o2));
		assertEquals(-1,comparator.compare(o2, o3));
		assertEquals(-1,comparator.compare(o1, o3));
		
		
		assertEquals(1,comparator.compare(o3, o2));
		assertEquals(1,comparator.compare(o2, o1));
		assertEquals(1,comparator.compare(o3, o1));
		
		assertEquals(0,comparator.compare(o1, o1));
		assertEquals(0,comparator.compare(o2, o2));
		assertEquals(0,comparator.compare(o3, o3));
		
		
		
	}
	

}
