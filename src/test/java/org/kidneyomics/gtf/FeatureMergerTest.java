package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;

public class FeatureMergerTest {

	@Test
	public void test() {
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(250, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f1);
		features.add(f2);
		features.add(f3);
		
		
		List<Feature> merged = FeatureMerger.mergeOverlappingFeatures(features);
		assertTrue(merged.size() == 2);
		
		Feature m1 = merged.get(0);
		
		assertEquals(50, m1.location().bioStart());
		assertEquals(200, m1.location().bioEnd());
		
		Feature m2 = merged.get(1);
		
		assertEquals(250, m2.location().bioStart());
		assertEquals(300, m2.location().bioEnd());
	}
	
	
	@Test
	public void test2() {
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '-'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(250, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f1);
		features.add(f2);
		features.add(f3);
		
		
		List<Feature> merged = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(features);
		assertTrue(merged.size() == 2);
		
		Feature m1 = merged.get(0);
		
		assertEquals(50, m1.location().bioStart());
		assertEquals(200, m1.location().bioEnd());
		
		Feature m2 = merged.get(1);
		
		assertEquals(250, m2.location().bioStart());
		assertEquals(300, m2.location().bioEnd());
	}
	
	@Test
	public void test3() {
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '-'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(250, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f1);
		features.add(f2);
		features.add(f3);
		
		IllegalArgumentException exception = null;
		try {
			FeatureMerger.mergeOverlappingFeatures(features);
		} catch(IllegalArgumentException e) {
			exception = e;
		}
		assertNotNull(exception);

	}

}
