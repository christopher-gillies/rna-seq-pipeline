package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;

public class FindOverlappingGenePairsTest {

	@Test
	public void test() {
		LinkedList<Feature> genes = new LinkedList<>();
		for(int i = 100; i <= 1000; i+=100) {
			Feature f = new Feature("chr1","","",Location.fromBio(i, i+200, '+'),0.0,0,"");
			genes.add(f);
		}
		/*
		 * [100,300]  ----- 	2
		 * 	[200,400] ----- 	2
		 *   [300,500] ---- 	2
		 *    [400,600]			2
		 *     [500,700]		2
		 *      [600,800]		2
		 *       [700,900]		2
		 *        [800, 1000]	2
		 *         [900, 1100]	1
		 *          [1000, 1200] 
		 */
		
		FindOverlappingGenePairs findOverlappingGenePairs = new FindOverlappingGenePairs(genes);
		List<FeaturePair> pairs = findOverlappingGenePairs.getOverlappingGenes();
		assertEquals(17,pairs.size());
		
		assertEquals(100, pairs.get(0).getFirst().location().bioStart());
		assertEquals(300, pairs.get(0).getFirst().location().bioEnd());
		
		assertEquals(200, pairs.get(0).getSecond().location().bioStart());
		assertEquals(400, pairs.get(0).getSecond().location().bioEnd());
		
		assertEquals(100, pairs.get(1).getFirst().location().bioStart());
		assertEquals(300, pairs.get(1).getFirst().location().bioEnd());
		
		assertEquals(300, pairs.get(1).getSecond().location().bioStart());
		assertEquals(500, pairs.get(1).getSecond().location().bioEnd());
		
		
		assertEquals(200, pairs.get(2).getFirst().location().bioStart());
		assertEquals(400, pairs.get(2).getFirst().location().bioEnd());
		
		assertEquals(300, pairs.get(2).getSecond().location().bioStart());
		assertEquals(500, pairs.get(2).getSecond().location().bioEnd());
		
		
		assertEquals(900, pairs.get(16).getFirst().location().bioStart());
		assertEquals(1100, pairs.get(16).getFirst().location().bioEnd());
		
		assertEquals(1000, pairs.get(16).getSecond().location().bioStart());
		assertEquals(1200, pairs.get(16).getSecond().location().bioEnd());
		
		for(FeaturePair pair : pairs) {
			assertTrue(pair.overlaps());
		}
	}

}
