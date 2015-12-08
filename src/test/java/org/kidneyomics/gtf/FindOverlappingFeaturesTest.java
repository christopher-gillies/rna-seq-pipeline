package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FindOverlappingFeaturesTest {

	
	Logger logger = LoggerFactory.getLogger(FindOverlappingFeaturesTest.class);
	@Test
	public void test() {
		
		
		Feature[] features = new Feature[100];
		
		int start = 1;
		int end = 50;
		for(int i = 0; i < 100; i++) {
			Location l1 = Location.fromBio(start, end, '+');
			
			start += 50;
			end += 50;
			
			features[i] = new Feature("chr1", "", "", l1, -1.0, -1, "");
			
			logger.info(GTFFeatureRenderer.render(features[i]));
			
		}
		
		Feature target = new Feature("chr1", "", "", Location.fromBio(400, 1000, '+'), -1.0, -1, "");
		
		int longest = FindOverlappingFeatures.findLongest(features);
		Arrays.sort(features, new FeatureComparator());
		List<Feature> overlap = FindOverlappingFeatures.findOverlappingFeatures(longest, features, target);

		logger.info("Overlapping features with " + GTFFeatureRenderer.render(target));
		for(Feature f : overlap) {
			logger.info(GTFFeatureRenderer.render(f));
		}
		
		assertEquals(13,overlap.size());
	
	}
	
	@Test
	public void test2() {
		
		logger.info("FindOverlappingFeaturesTest#test2");
		
		Feature[] features = new Feature[100];
		
		int start = 1;
		int end = 50;
		for(int i = 0; i < 100; i++) {
			Location l1 = Location.fromBio(start, end, '+');
			
			start += 50;
			end += 50 + i;
			
			features[i] = new Feature("chr1", "", "", l1, -1.0, -1, "");
			
			logger.info(GTFFeatureRenderer.render(features[i]));
			
		}
		
		Feature target = new Feature("chr1", "", "", Location.fromBio(350, 1001, '+'), -1.0, -1, "");
		
		int longest = FindOverlappingFeatures.findLongest(features);
		Arrays.sort(features, new FeatureComparator());
		List<Feature> overlap = FindOverlappingFeatures.findOverlappingFeatures(longest, features, target);

		logger.info("Overlapping features with " + GTFFeatureRenderer.render(target));
		for(Feature f : overlap) {
			logger.info(GTFFeatureRenderer.render(f));
		}
		
		assertEquals(15,overlap.size());
	
	}
	
	
	@Test
	public void test3() {
		
		logger.info("FindOverlappingFeaturesTest#test3");
		
		Feature[] features = new Feature[100];
		
		int start = 1;
		int end = 50;
		for(int i = 0; i < 100; i++) {
			Location l1 = Location.fromBio(start, end, '+');
			
			start += 50;
			end += 50 + i;
			
			features[i] = new Feature("chr1", "", "", l1, -1.0, -1, "");
			
			logger.info(GTFFeatureRenderer.render(features[i]));
			
		}
		
		Feature target = new Feature("chr1", "", "", Location.fromBio(925, 1001, '+'), -1.0, -1, "");
		
		int longest = FindOverlappingFeatures.findLongest(features);
		Arrays.sort(features, new FeatureComparator());
		List<Feature> overlap = FindOverlappingFeatures.findOverlappingFeatures(longest, features, target);

		logger.info("Overlapping features with " + GTFFeatureRenderer.render(target));
		for(Feature f : overlap) {
			logger.info(GTFFeatureRenderer.render(f));
		}
		
		assertEquals(5,overlap.size());
	
	}

}
