package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FeatureMergerTest {

	
	Logger logger = LoggerFactory.getLogger(FeatureMergerTest.class);
	
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
	
	@Test
	public void test4() {
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(1, 400, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '-'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(250, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f1);
		features.add(f2);
		features.add(f3);
		
		
		List<Feature> merged = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(features);
		assertTrue(merged.size() == 1);
		
		Feature m1 = merged.get(0);
		
		assertEquals(1, m1.location().bioStart());
		assertEquals(400, m1.location().bioEnd());
		
	}
	
	
	@Test
	public void testRemoveOverlappingFeatures() {
		
		Feature f0 = new Feature("chr1", "a", "exon", Location.fromBio(50, 99, '+'), 0.0, 0, "");
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(150, 250, '+'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(260, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f0);
		features.add(f1);
		features.add(f2);
		features.add(f3);
		
		List<Feature> removed = FeatureMerger.removeOverlappingFeaturesIgnoreStrand(features);
		
		logger.info("removed");
		for(Feature f : removed) {
			logger.info(GTFFeatureRenderer.render(f));
		}
		
		logger.info("kept");
		for(Feature f : features) {
			logger.info(GTFFeatureRenderer.render(f));
		}
		
		assertEquals(2,removed.size());
		assertEquals(2,features.size());
		
		assertEquals(100,removed.get(0).location().bioStart());
		assertEquals(200,removed.get(0).location().bioEnd());
		
		assertEquals(150,removed.get(1).location().bioStart());
		assertEquals(250,removed.get(1).location().bioEnd());
		
		assertEquals(50,features.get(0).location().bioStart());
		assertEquals(99,features.get(0).location().bioEnd());
		
		assertEquals(260,features.get(1).location().bioStart());
		assertEquals(300,features.get(1).location().bioEnd());
		
	}
	
	@Test
	public void testRemoveOverlappingFeatures2() {
		
		Feature f0 = new Feature("chr1", "a", "exon", Location.fromBio(1, 261, '+'), 0.0, 0, "");
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(150, 250, '+'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(260, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f0);
		features.add(f1);
		features.add(f2);
		features.add(f3);
		
		List<Feature> removed = FeatureMerger.removeOverlappingFeaturesIgnoreStrand(features);
		
		assertEquals(4,removed.size());
		assertEquals(0,features.size());
		

		
		assertEquals(100,removed.get(0).location().bioStart());
		assertEquals(200,removed.get(0).location().bioEnd());
		
		assertEquals(150,removed.get(1).location().bioStart());
		assertEquals(250,removed.get(1).location().bioEnd());
		
		assertEquals(1,removed.get(2).location().bioStart());
		assertEquals(261,removed.get(2).location().bioEnd());
		
		assertEquals(260,removed.get(3).location().bioStart());
		assertEquals(300,removed.get(3).location().bioEnd());
		
	}
	
	@Test
	public void testRemoveOverlappingFeatures3() {
		
		Feature f0 = new Feature("chr1", "a", "exon", Location.fromBio(1, 50, '+'), 0.0, 0, "");
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 148, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(150, 250, '+'), 0.0, 0, "");
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(260, 400, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr1", "a", "exon", Location.fromBio(260, 300, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f0);
		features.add(f1);
		features.add(f2);
		features.add(f3);
		features.add(f4);
		
		List<Feature> removed = FeatureMerger.removeOverlappingFeaturesIgnoreStrand(features);
		
		assertEquals(2,removed.size());
		assertEquals(3,features.size());
		

		
		assertEquals(260,removed.get(0).location().bioStart());
		assertEquals(400,removed.get(0).location().bioEnd());
		
		assertEquals(260,removed.get(1).location().bioStart());
		assertEquals(300,removed.get(1).location().bioEnd());
		
		assertEquals(1,features.get(0).location().bioStart());
		assertEquals(50,features.get(0).location().bioEnd());
		
		assertEquals(100,features.get(1).location().bioStart());
		assertEquals(148,features.get(1).location().bioEnd());
		
		assertEquals(150,features.get(2).location().bioStart());
		assertEquals(250,features.get(2).location().bioEnd());
		
	}

	@Test
	public void testRemoveOverlappingFeatures4() {
		
		Feature f0 = new Feature("chr1", "a", "exon", Location.fromBio(1, 261, '+'), 0.0, 0, "");
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(150, 250, '+'), 0.0, 0, "");
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(260, 300, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr1", "a", "exon", Location.fromBio(300, 400, '+'), 0.0, 0, "");
		
		Feature f5 = new Feature("chr1", "a", "exon", Location.fromBio(500, 600, '+'), 0.0, 0, "");
		
		List<Feature> features = new LinkedList<Feature>();
		
		features.add(f0);
		features.add(f1);
		features.add(f2);
		features.add(f3);
		features.add(f4);
		features.add(f5);
		
		List<Feature> removed = FeatureMerger.removeOverlappingFeaturesIgnoreStrand(features);
		
		assertEquals(5,removed.size());
		assertEquals(1,features.size());
		

		
		assertEquals(100,removed.get(0).location().bioStart());
		assertEquals(200,removed.get(0).location().bioEnd());
		
		assertEquals(150,removed.get(1).location().bioStart());
		assertEquals(250,removed.get(1).location().bioEnd());
		
		assertEquals(1,removed.get(2).location().bioStart());
		assertEquals(261,removed.get(2).location().bioEnd());
		
		assertEquals(260,removed.get(3).location().bioStart());
		assertEquals(300,removed.get(3).location().bioEnd());
		
		assertEquals(300,removed.get(4).location().bioStart());
		assertEquals(400,removed.get(4).location().bioEnd());
		
		assertEquals(500,features.get(0).location().bioStart());
		assertEquals(600,features.get(0).location().bioEnd());
		
	}
}