package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class LocationComparatorTest {

	
	Logger logger = LoggerFactory.getLogger(LocationComparatorTest.class);
	
	@Test
	public void test() {
		LocationComparator comparator = new LocationComparator();
		
		Location o1 = Location.fromBio(100, 200, '+');
		Location o2 = Location.fromBio(150, 200, '-');
		Location o3 = Location.fromBio(160, 200, '-');
		
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
