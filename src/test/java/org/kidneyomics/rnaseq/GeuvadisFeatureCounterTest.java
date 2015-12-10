package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.kidneyomics.gtf.FeatureMerger;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.SAMRecordToFeatureConverter;
import org.mockito.Matchers;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import static org.mockito.Mockito.*;

public class GeuvadisFeatureCounterTest {

	Logger logger = LoggerFactory.getLogger(GeuvadisFeatureCounterTest.class);
	
	@Test
	public void testBuildFeatures() throws IOException {
		
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		
		File gtf = r.getFile();
		
		HashMap<String,Integer> genes = new HashMap<>();
		
		HashSet<Integer> totalCoverage = new HashSet<>();
		
		int longestExon = 0;
		try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
			for(Feature f : reader) {
				if(f.type().equals("exon")) {
					String geneId = f.getAttribute("gene_id");
					genes.put(geneId, 0);
					
					if(f.location().length() > longestExon) {
						logger.info("Feature:\n" + GTFFeatureRenderer.render(f));
						longestExon = f.location().length();
					}
					
					for(int i = f.location().bioStart(); i <= f.location().bioEnd(); i++) {
						totalCoverage.add(i);
					}
				}
			}
		}
		
		logger.info("Longest feature: " + longestExon);
		
		logger.info("Total coverage: " + totalCoverage.size());
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GeuvadisFeatureCounter gfc = new GeuvadisFeatureCounter(findOverlappingFeatures, new LoggerService());
		
		gfc.buildFeatures(gtf, "exon");
		
		
		logger.info( "Number of chromosomes: " + gfc.getChromosomeFeatures().size() );
		
		logger.info("Number of merged features per gene across chromosome 1: " + gfc.getChromosomeFeatures().get("chr1").length);
		
		logger.info("Longest feature on chromosme 1: " + gfc.getChromosomeLongestFeature().get("chr1"));
		
		//validate that the features are sorted
		Feature[] features = gfc.getChromosomeFeatures().get("chr1");
		for(int i = 0; i < features.length - 1; i++) {
			assertTrue(features[i].location().bioStart() <= features[i+1].location().bioStart());
		}
		
		
		for(Feature f : features) {
			String geneId = f.getAttribute("gene_id");
			assertTrue(genes.containsKey(geneId));
			
			genes.put(geneId, 1);
		}
		
		//validate all genes are present
		for(Map.Entry<String, Integer> entry : genes.entrySet()) {
			assertTrue(entry.getValue() > 0);
		}
		
		//make sure the same total coverage is found
		int lengthDoubleCountingOverlappingGenes = 0;
		int longestGeneLevelMerged = 0;
		for(int i = 0; i < features.length; i++) {
			lengthDoubleCountingOverlappingGenes += features[i].location().length();
			
			if(features[i].location().length() > longestGeneLevelMerged) {
				longestGeneLevelMerged = features[i].location().length();
			}
		}
		
		logger.info("Longest feature merged: " + longestGeneLevelMerged);
		
		assertTrue(longestGeneLevelMerged >= longestExon);
		
		logger.info("Total coverage calculated double counting overlapping genes: " + lengthDoubleCountingOverlappingGenes);
		
		assertTrue(lengthDoubleCountingOverlappingGenes >= totalCoverage.size());
		
		
		List<Feature> mergedFeatures = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(Arrays.asList(features));
		int length = 0;
		for(Feature f : mergedFeatures) {
			length += f.location().length();
		}
		logger.info("Total coverage: " + length);
		
		assertEquals(totalCoverage.size(),length);
		
	}

	@Test
	public void testmapToMultipleGenes() {
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GeuvadisFeatureCounter gfc = new GeuvadisFeatureCounter(findOverlappingFeatures, new LoggerService());
		String line1 = "chr1	HAVANA	exon	788771	794826	.	+	.	gene_id \"ENSG00000228794\"; transcript_id \"ENST00000445118\";";
		Feature f1 = GTFFeatureBuilder.createFromLine(line1);
		
		String line2 = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG00000228794\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		
		HashSet<Feature> features = new HashSet<>();
		
		features.add(f1);
		features.add(f1);
		
		assertTrue(features.size() == 1);
		
		features.add(f2);
		assertTrue(features.size() == 2);
		
		assertTrue(gfc.mapToMultipleGenes(features) == false);
		
		
		String line3 = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f3 = GTFFeatureBuilder.createFromLine(line3);
		
		features.add(f3);
		assertTrue(features.size() == 3);
		
		assertTrue(gfc.mapToMultipleGenes(features) == true);
	}
	
	@Test
	public void test1GetMappedRegionsForMate() {
		
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		record.setCigarString("39M");
		
		SAMRecordToFeatureConverter converter = new SAMRecordToFeatureConverter();
		List<Feature> featuresForRecord = converter.convert(record);
		assertTrue(featuresForRecord.size() == 1);
		assertTrue(featuresForRecord.get(0).location().length() == 39);
		assertEquals("chr1",featuresForRecord.get(0).seqname());
		assertEquals(100,featuresForRecord.get(0).location().bioStart());
		assertEquals(138,featuresForRecord.get(0).location().bioEnd());
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		
		LinkedList<Feature> result = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result.add(f);
		
		when(mock.findOverlappingFeatures(anyInt(), any(Feature[].class), any(Feature.class))).thenReturn(result);
		
		GeuvadisFeatureCounter gfc = new GeuvadisFeatureCounter(mock, new LoggerService());
		
		Set<Feature> featuresForMate = gfc.getMappedRegionsForMate(record, -1, null);
		
		
		assertEquals(1,featuresForMate.size());
	}
	
	@Test
	public void test2GetMappedRegionsForMate() {
		
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		record.setCigarString("39M");
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		

		LinkedList<Feature> result = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result.add(f);
		
		String line2 = "chr1	HAVANA	exon	788771	794826	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		result.add(f2);
		
		when(mock.findOverlappingFeatures(anyInt(), any(Feature[].class), any(Feature.class))).thenReturn(result);
		
		GeuvadisFeatureCounter gfc = new GeuvadisFeatureCounter(mock, new LoggerService());
		
		Set<Feature> featuresForMate = gfc.getMappedRegionsForMate(record, -1, null);
				
		assertEquals(2,featuresForMate.size());
	}
	
	@Test
	public void test3GetMappedRegionsForMate() {
		
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		record.setCigarString("39M10N1M");
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		

		LinkedList<Feature> result1 = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result1.add(f);
		
		String line2 = "chr1	HAVANA	exon	788771	794826	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		result1.add(f2);
		
		LinkedList<Feature> result2 = new LinkedList<Feature>();
		
		String line3 = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f3 = GTFFeatureBuilder.createFromLine(line3);
		result2.add(f3);
		
		
		when(mock.findOverlappingFeatures(anyInt(), any(Feature[].class), any(Feature.class))).thenReturn(result1,result2);
		
		GeuvadisFeatureCounter gfc = new GeuvadisFeatureCounter(mock, new LoggerService());
		
		Set<Feature> featuresForMate = gfc.getMappedRegionsForMate(record, -1, null);
				
		assertEquals(3,featuresForMate.size());
	}
	
	
	@Test
	public void test4GetMappedRegionsForMate() {
		
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		record.setCigarString("39M10N1M");
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		

		LinkedList<Feature> result1 = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result1.add(f);
		
		String line2 = "chr1	HAVANA	exon	788771	794826	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		result1.add(f2);
		
		LinkedList<Feature> result2 = new LinkedList<Feature>();
		
		String line3 = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f3 = GTFFeatureBuilder.createFromLine(line3);
		result2.add(f3);
		
		
		when(mock.findOverlappingFeatures(anyInt(), any(Feature[].class), any(Feature.class))).thenReturn(result1,result2,result1);
		
		GeuvadisFeatureCounter gfc = new GeuvadisFeatureCounter(mock, new LoggerService());
		
		Set<Feature> featuresForMate = gfc.getMappedRegionsForMate(record, -1, null);
				
		assertEquals(3,featuresForMate.size());
	}
	
	@Test
	public void testCountUnpaired() {
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		record.setCigarString("39M10N1M");
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		
		
		/*
		 * test counts
		 */
		when(mock.findOverlappingFeatures(anyInt(), any(Feature[].class), any(Feature.class))).thenReturn(null);
		
	}
}
