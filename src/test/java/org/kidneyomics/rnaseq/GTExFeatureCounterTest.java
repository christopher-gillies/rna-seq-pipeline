package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.kidneyomics.gtf.ExonFilter;
import org.kidneyomics.gtf.FeatureComparator;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.FeatureMerger;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFFeatureUtil;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.RemoveRetainedIntronFilter;
import org.kidneyomics.gtf.SAMRecordToFeatureConverter;
import org.mockito.Matchers;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.Assert;

import static org.mockito.Mockito.*;

public class GTExFeatureCounterTest {

	Logger logger = LoggerFactory.getLogger(GTExFeatureCounterTest.class);
	
	
	
	@Test
	public void testBuildFeatures() throws IOException {
		
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		
		List<Feature> allExons = new LinkedList<>();
		
		File gtf = r.getFile();
		
		HashMap<String,Integer> genes = new HashMap<>();
		
		
		try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
			//Only read exons
			reader.addFilter(new ExonFilter()).addFilter(new RemoveRetainedIntronFilter());
			for(Feature f : reader) {
				String geneId = f.getAttribute("gene_id");
				genes.put(geneId, 0);
				allExons.add(f);
			}
		}
		
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		
		gfc.buildFeatures(gtf, "exon");
		
		HashSet<Feature> removedFeatures = gfc.getRemovedFeatures();
		
		logger.info( "Number of chromosomes: " + gfc.getChromosomeFeatures().size() );
		
		logger.info("Number of merged features per gene across chromosome 1: " + gfc.getChromosomeFeatures().get("chr1").length);
		
		//validate that the features are sorted
		Feature[] features = gfc.getChromosomeFeatures().get("chr1");
		assertTrue(GTFFeatureUtil.isSorted(features));
		assertTrue(GTFFeatureUtil.hasNoOverlapIgnoreStrand(features));
		
		logger.info("features sorted");
		
		
		for(Feature f : features) {
			String geneId = f.getAttribute("gene_id");
			assertTrue(genes.containsKey(geneId));
			
			genes.put(geneId, 1);
		}
		
		
		
		logger.info(removedFeatures.size() + " Remove features ");
		HashMap<String,Feature> removedGeneNames = new HashMap<>();

		
		logger.warn("GENES NOT FOUND:");
		for(Feature f : removedFeatures) {
			removedGeneNames.put(f.getAttribute("gene_id"),f);
			logger.info(f.getAttribute("gene_id"));
			//logger.info(GTFFeatureRenderer.render(f));
		}
		
		
		//validate all genes are present
		int removedGenes = 0;
		int totalGenes = 0;
		for(Map.Entry<String, Integer> entry : genes.entrySet()) {
			totalGenes++;
			if(entry.getValue() == 0) {
				//logger.warn(entry.getKey());
				removedGenes++;
				assertTrue(removedGeneNames.containsKey(entry.getKey()));
				
			}
		}
		logger.warn("Removed genes: " + removedGenes + " of " + totalGenes);
		
		//make sure the same total coverage is found
		int totalCoveragePreSecondMerge = 0;
		int longestGeneLevelMerged = 0;
		for(int i = 0; i < features.length; i++) {
			totalCoveragePreSecondMerge += features[i].location().length();
			
			if(features[i].location().length() > longestGeneLevelMerged) {
				longestGeneLevelMerged = features[i].location().length();
			}
		}
		logger.info("Total coverage before merge: " + totalCoveragePreSecondMerge);
		
		
		
		List<Feature> mergedFeatures = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(Arrays.asList(features));
		int length = 0;
		for(Feature f : mergedFeatures) {
			length += f.location().length();
		}
		logger.info("Total coverage: " + length);
		
		//assertEquals(totalCoveragePreSecondMerge,length);
		
		
		HashSet<Integer> totalCoverage = new HashSet<>();
		for(Feature f : features) {
			if(!removedFeatures.contains(f)) {
				for(int i = f.location().bioStart(); i <= f.location().bioEnd(); i++) {
					totalCoverage.add(i);
				}
			}
		}
		logger.info("Total coverage calculated excluding multigene features: " + totalCoverage.size());
		
		//validate feature counts
		Collection<String> fcs = gfc.getFeatureIds();
		
		assertEquals(gfc.getChromosomeFeatures().get("chr1").length, fcs.size());
		
		for(String id : fcs) {
			FeatureCount fc = gfc.getCounts(id);
			String fcId = fc.getId();
			String featureId = fc.getFeature().getAttribute("id");
			
			assertEquals(id,fcId);
			assertEquals(id,featureId);
			//
		}
		
		List<FeatureCount> fcList = gfc.getCounts();
		Collections.sort(fcList);
		for(FeatureCount fc : fcList) {
			//logger.info(fc.getId());
		}
	}

	@Test
	public void testMapToMultipleGenes() {
		logger.info("testMapToMultipleGenes");
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
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
		logger.info("test1GetMappedRegionsForMate");
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(788771);
		record.setCigarString("39M");
		
		SAMRecordToFeatureConverter converter = new SAMRecordToFeatureConverter();
		List<Feature> featuresForRecord = converter.convert(record);
		assertTrue(featuresForRecord.size() == 1);
		assertTrue(featuresForRecord.get(0).location().length() == 39);
		assertEquals("chr1",featuresForRecord.get(0).seqname());
		assertEquals(788771,featuresForRecord.get(0).location().bioStart());
		assertEquals(788771 + 38,featuresForRecord.get(0).location().bioEnd());
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		
		LinkedList<Feature> result = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	788771	794827	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result.add(f);
		
		when(mock.findOverlappingFeatures(any(Feature[].class), any(Feature.class))).thenReturn(result);
		
		GTExFeatureCounter gfc = new GTExFeatureCounter(mock, new LoggerService());
		
		Map<Feature,Integer> featuresForMate = gfc.getMappedRegionsForMate(record,  null);
		
		
		assertEquals(1,featuresForMate.size());
		assertEquals(39,featuresForMate.get(f).intValue());
	}
	
	@Test
	public void test2GetMappedRegionsForMate() {
		logger.info("test2GetMappedRegionsForMate");
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(150);
		record.setCigarString("51M100N30M");
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		

		LinkedList<Feature> result = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	100	200	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result.add(f);
		
		String line2 = "chr1	HAVANA	exon	300	400	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		result.add(f2);
		
		when(mock.findOverlappingFeatures(any(Feature[].class), any(Feature.class))).thenReturn(result);
		
		GTExFeatureCounter gfc = new GTExFeatureCounter(mock, new LoggerService());
		
		Map<Feature,Integer> featuresForMate = gfc.getMappedRegionsForMate(record,  null);
		
		
		assertEquals(2,featuresForMate.size());
		assertEquals(51,featuresForMate.get(f).intValue());
		assertEquals(30,featuresForMate.get(f2).intValue());
	}
	
	//Test some percentage not mapping to the exon
	@Test
	public void test3GetMappedRegionsForMate() {
		logger.info("test3GetMappedRegionsForMate");
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(100);
		record.setCigarString("10M9N10M9N10M");
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		

		LinkedList<Feature> result1 = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	100	110	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result1.add(f);
		
		String line2 = "chr1	HAVANA	exon	120	130	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		result1.add(f2);
		
		String line3 = "chr1	HAVANA	exon	140	150	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f3 = GTFFeatureBuilder.createFromLine(line3);
		result1.add(f3);
		
		
		when(mock.findOverlappingFeatures( any(Feature[].class), any(Feature.class))).thenReturn(result1);
		
		GTExFeatureCounter gfc = new GTExFeatureCounter(mock, new LoggerService());
		
		Map<Feature,Integer> featuresForMate = gfc.getMappedRegionsForMate(record,  null);
		
		
		assertEquals(3,featuresForMate.size());
		assertEquals(10,featuresForMate.get(f).intValue());
		assertEquals(9,featuresForMate.get(f2).intValue());
		assertEquals(8,featuresForMate.get(f3).intValue());
	}
	
	@Test
	public void testAddToFeatures() {
		logger.info("testAddToFeatures");
		
		//static void addToFeatures(Map<Feature,Integer> features, Map<String,FeatureCount> featureCounts) {
		
		
		
		String line = "chr1	HAVANA	exon	100	110	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f1 = GTFFeatureBuilder.createFromLine(line);
		FeatureCount fc1 = new FeatureCount(f1);
		
		String line2 = "chr1	HAVANA	exon	120	130	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		FeatureCount fc2 = new FeatureCount(f2);
		
		String line3 = "chr1	HAVANA	exon	140	150	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f3 = GTFFeatureBuilder.createFromLine(line3);
		FeatureCount fc3 = new FeatureCount(f3);
		
		HashMap<Feature,Integer> features = new HashMap<>();
		features.put(f1, 10);
		features.put(f2, 11);
		features.put(f3, 12);
		
		HashMap<String,FeatureCount> fcCounts = new HashMap<>();
		fcCounts.put(fc1.getId(), fc1);
		fcCounts.put(fc2.getId(), fc2);
		fcCounts.put(fc3.getId(), fc3);
		
		GTExFeatureCounter.addToFeatures(features, fcCounts);
		
		assertEquals(10.0 / (33.0),  fcCounts.get(fc1.getId()).getCount(), 0.0001);
		assertEquals(11.0 / (33.0),  fcCounts.get(fc2.getId()).getCount(), 0.0001);
		assertEquals(12.0 / (33.0),  fcCounts.get(fc3.getId()).getCount(), 0.0001);
	}
	
	@Test
	public void testCountPaire1() throws IOException {
		logger.info("testCountPaire1");
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		record1.setReferenceName("chr1");
		record1.setAlignmentStart(11869);
		record1.setCigarString("39M");
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(12000);
		record2.setCigarString("39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		
		assertEquals(1.0,fc.getCount(),0.000001);
		
		assertEquals(1,gfc.getTotalCount());
		assertEquals(1,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(0,gfc.getAmbiguousReadCount());
	}
	
	
	//@Test
	public void testCountPaire2() throws IOException {
		logger.info("testCountPaire1");
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		record1.setReferenceName("chr1");
		record1.setAlignmentStart(11869);
		record1.setCigarString("39M");
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(12000);
		record2.setCigarString("39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		
		assertEquals(1.0,fc.getCount(),0.000001);
		
		assertEquals(1,gfc.getTotalCount());
		assertEquals(1,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(0,gfc.getAmbiguousReadCount());
	}
	
	
	//@Test
	public void testCountUnpaired2() throws IOException {
		
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(11869);
		record.setCigarString("39M800N39M250N39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record);
		
		gfc.count(pair);
		
		assertEquals(1/3.0,gfc.getCounts("ENSG00000223972_chr1_11869_12227").getCount(),0.000001);
		assertEquals(1/3.0,gfc.getCounts("ENSG00000223972_chr1_12595_12721").getCount(),0.000001);
		assertEquals(1/3.0,gfc.getCounts("ENSG00000223972_chr1_12975_13052").getCount(),0.000001);
		
		assertEquals(1,gfc.getTotalCount());
		assertEquals(1,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(0,gfc.getAmbiguousReadCount());
		
		gfc.count(pair);
		
		assertEquals(2/3.0,gfc.getCounts("ENSG00000223972_chr1_11869_12227").getCount(),0.000001);
		assertEquals(2/3.0,gfc.getCounts("ENSG00000223972_chr1_12595_12721").getCount(),0.000001);
		assertEquals(2/3.0,gfc.getCounts("ENSG00000223972_chr1_12975_13052").getCount(),0.000001);
		
		assertEquals(2,gfc.getTotalCount());
		assertEquals(2,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(0,gfc.getAmbiguousReadCount());
		
		gfc.count(pair);
		
		assertEquals(1.0,gfc.getCounts("ENSG00000223972_chr1_11869_12227").getCount(),0.000001);
		assertEquals(1.0,gfc.getCounts("ENSG00000223972_chr1_12595_12721").getCount(),0.000001);
		assertEquals(1.0,gfc.getCounts("ENSG00000223972_chr1_12975_13052").getCount(),0.000001);
		
		assertEquals(3,gfc.getTotalCount());
		assertEquals(3,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(0,gfc.getAmbiguousReadCount());
		
	}
	
	//@Test
	public void testCountAmbiguous() throws IOException {
		
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(11869);
		record.setCigarString("39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		
		assertEquals(1.0,fc.getCount(),0.000001);
		
		assertEquals(1,gfc.getTotalCount());
		assertEquals(1,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(0,gfc.getAmbiguousReadCount());
		
		
		SAMRecordPair pair2 = new SAMRecordPair();
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(11869);
		record2.setCigarString("39M17685N39M");
		
		pair2.setMate1(record2);
		
		gfc.count(pair2);
		
		assertEquals(2,gfc.getTotalCount());
		assertEquals(1,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(1,gfc.getAmbiguousReadCount());
		
		// Second ambiguous
		
		//ENSG00000238009_chr1_89295_91629
		//ENSG00000239945_chr1_89551_90050
		
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(89551);
		record2.setCigarString("39M");
		gfc.count(pair2);
		
		assertEquals(3,gfc.getTotalCount());
		assertEquals(1,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(2,gfc.getAmbiguousReadCount());
	}
	
	
	//@Test
	public void testCountAmbiguousPairOnePairAmbiguous() throws IOException {
		
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(11869);
		record.setCigarString("39M");
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(11869);
		record2.setCigarString("39M17685N39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record);
		pair.setMate2(record2);
		
		
		gfc.count(pair);
		
		
		assertEquals(2,gfc.getTotalCount());
		assertEquals(0,gfc.getMappedReadCount());
		assertEquals(0,gfc.getUnmappedReadCount());
		assertEquals(2,gfc.getAmbiguousReadCount());
		
	}
}
