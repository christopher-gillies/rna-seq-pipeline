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
		
		logger.info("FindOverlappingFeaturesTest#test");
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
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		
		Arrays.sort(features, new FeatureComparator());
		List<Feature> overlap = findOverlappingFeatures.findOverlappingFeatures(features, target);

		logger.info("Overlapping features with " + GTFFeatureRenderer.render(target));
		for(Feature f : overlap) {
			logger.info(GTFFeatureRenderer.render(f));
		}
		
		assertEquals(13,overlap.size());
	
	}
	
}
