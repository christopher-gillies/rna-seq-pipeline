package org.kidneyomics.rnaseq;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.Matchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.kidneyomics.gtf.ExonFilter;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.FeatureMerger;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFFeatureUtil;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.RemoveRetainedIntronFilter;
import org.kidneyomics.gtf.SAMRecordToFeatureConverter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

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
		
		//List<FeatureCount> fcList = gfc.getCounts();
		//Collections.sort(fcList);
		//for(FeatureCount fc : fcList) {
			//logger.info(fc.getId());
		//}
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
		
		Map<Feature,Integer> featuresForMate = GTExFeatureCounter.getMappedRegionsForMate(record,  null, new SAMRecordToFeatureConverter(), mock);
		
		
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
		
		Map<Feature,Integer> featuresForMate = GTExFeatureCounter.getMappedRegionsForMate(record,  null, new SAMRecordToFeatureConverter(), mock);
		
		
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
		
		Map<Feature,Integer> featuresForMate = GTExFeatureCounter.getMappedRegionsForMate(record,  null, new SAMRecordToFeatureConverter(), mock);
		
		
		assertEquals(3,featuresForMate.size());
		assertEquals(10,featuresForMate.get(f).intValue());
		assertEquals(9,featuresForMate.get(f2).intValue());
		assertEquals(8,featuresForMate.get(f3).intValue());
	}
	
	
	@Test
	public void test4GetMappedRegionsForMate() {
		logger.info("test4GetMappedRegionsForMate");
		//getMappedRegionsForMate(SAMRecord mate, Set<Feature> featuresForMate, int longest, Feature[] features)
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		//D7DHSVN1:224:C28PGACXX:7:1302:20524:70946	1107	1	565317	255	14M1D25M	=	565067	-290	ATGGCTATAGCAATAAACTAGGAATAGCCCCCTTTCACT	=8F=9/*894?D9<JIIJJJJJJJJJJJJJJJJJIHFJH	PG:Z:MarkDuplicates	RG:Z:25979	NH:i:1	HI:i:1	nM:i:0	AS:i:72

		record.setReferenceName("chr1");
		record.setAlignmentStart(565317);
		record.setCigarString("14M1D25M");
		
		SAMRecordToFeatureConverter converter = new SAMRecordToFeatureConverter();
		List<Feature> featuresForRecord = converter.convert(record);
		assertTrue(featuresForRecord.size() == 2);
		assertTrue(featuresForRecord.get(0).location().length() == 14);
		assertTrue(featuresForRecord.get(1).location().length() == 25);
		assertEquals("chr1",featuresForRecord.get(0).seqname());
		assertEquals(565317,featuresForRecord.get(0).location().bioStart());
		assertEquals(565317 + 13,featuresForRecord.get(0).location().bioEnd());
		
		
		assertEquals("chr1",featuresForRecord.get(1).seqname());
		assertEquals(565317 + 13 + 2,featuresForRecord.get(1).location().bioStart());
		assertEquals(565317 + 13 + 2 + 24,featuresForRecord.get(1).location().bioEnd());
		
		FindOverlappingFeatures mock = mock(FindOverlappingFeatures.class);
		
		LinkedList<Feature> result = new LinkedList<Feature>();
		
		String line = "chr1	HAVANA	exon	565317	566317	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f = GTFFeatureBuilder.createFromLine(line);
		result.add(f);
		
		when(mock.findOverlappingFeatures(any(Feature[].class), any(Feature.class))).thenReturn(result);
		
		Map<Feature,Integer> featuresForMate = GTExFeatureCounter.getMappedRegionsForMate(record,  null, new SAMRecordToFeatureConverter(), mock);
		
		
		assertEquals(1,featuresForMate.size());
		assertEquals(39,featuresForMate.get(f).intValue());
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
		
		GTExFeatureCounter.addToFeatures(100, features, fcCounts);
		
		assertEquals(10.0 / (100.0),  fcCounts.get(fc1.getId()).getCount(), 0.0001);
		assertEquals(11.0 / (100.0),  fcCounts.get(fc2.getId()).getCount(), 0.0001);
		assertEquals(12.0 / (100.0),  fcCounts.get(fc3.getId()).getCount(), 0.0001);
	}
	
	@Test
	public void testUnionFeaturesForMates() {
		logger.info("testUnionFeaturesForMates");
		//static Map<Feature,Integer> unionFeaturesForMates(Map<Feature,Integer> mappedFeaturesAcrossChunksForMate1, Map<Feature,Integer> mappedFeaturesAcrossChunksForMate2)
		
		String line = "chr1	HAVANA	exon	100	110	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f1 = GTFFeatureBuilder.createFromLine(line);
		
		String line2 = "chr1	HAVANA	exon	120	130	.	+	.	gene_id \"ENSG0000022879\"; transcript_id \"ENST00000445118\";";
		Feature f2 = GTFFeatureBuilder.createFromLine(line2);
		
		Map<Feature,Integer> mate1 = new HashMap<>();
		mate1.put(f1, 3);
		mate1.put(f2, 5);
		
		
		Map<Feature,Integer> mate2 = new HashMap<>();
		mate2.put(f1, 7);
		mate2.put(f2, 10);
		
		Map<Feature,Integer> union = GTExFeatureCounter.unionFeaturesForMates(mate1, mate2);
		
		assertEquals(2,union.size());
		assertEquals(10,union.get(f1).intValue());
		assertEquals(15,union.get(f2).intValue());
		
	}
	
	@Test
	public void testCountPair1() throws IOException {
		logger.info("testCountPair1");
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
		logger.info(GTFFeatureRenderer.render(fc.getFeature()));
		
		logger.info("Total counts for ENSG00000223972_chr1_11869_12227: " + fc.getCount());
		logger.info("Total counts reads: " + gfc.getTotalCount());
		logger.info("Total counts mappedReadCount: " + gfc.getMappedReadCount());
		logger.info("Total counts unmappedReadCount: " + gfc.getUnmappedReadCount());
		logger.info("Total counts ambiguous read count: " + gfc.getAmbiguousReadCount());
		logger.info("Partially unmapped reads: " + gfc.getNumberOfPartiallyUnmappedReads());
		
		assertEquals(1.0,fc.getCount(),0.000001);
		
		assertEquals(1.0,gfc.getTotalCount(),0.000001);
		assertEquals(1.0,gfc.getMappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(0,gfc.getNumberOfPartiallyUnmappedReads());
	}
	
	
	@Test
	public void testCountPair2() throws IOException {
		logger.info("testCountPair1");
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		record1.setReferenceName("chr1");
		record1.setAlignmentStart(11869);
		record1.setCigarString("50M");
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(12200);
		record2.setCigarString("50M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		logger.info(GTFFeatureRenderer.render(fc.getFeature()));
		
		logger.info("Total counts for ENSG00000223972_chr1_11869_12227: " + fc.getCount());
		logger.info("Total counts reads: " + gfc.getTotalCount());
		logger.info("Total counts mappedReadCount: " + gfc.getMappedReadCount());
		logger.info("Total counts unmappedReadCount: " + gfc.getUnmappedReadCount());
		logger.info("Total counts ambiguous read count: " + gfc.getAmbiguousReadCount());
		logger.info("Partially unmapped reads: " + gfc.getNumberOfPartiallyUnmappedReads());
		
		assertEquals( 78 / 100.0,fc.getCount(),0.000001);
		
		assertEquals(1.0,gfc.getTotalCount(),0.000001);
		assertEquals(78 / 100.0,gfc.getMappedReadCount(),0.000001);
		assertEquals(22 / 100.0,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(1,gfc.getNumberOfPartiallyUnmappedReads());
		
		
		/*
		 * Count again
		 */
		
		gfc.count(pair);
		
		fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		logger.info(GTFFeatureRenderer.render(fc.getFeature()));
		
		logger.info("Total counts for ENSG00000223972_chr1_11869_12227: " + fc.getCount());
		logger.info("Total counts reads: " + gfc.getTotalCount());
		logger.info("Total counts mappedReadCount: " + gfc.getMappedReadCount());
		logger.info("Total counts unmappedReadCount: " + gfc.getUnmappedReadCount());
		logger.info("Total counts ambiguous read count: " + gfc.getAmbiguousReadCount());
		logger.info("Partially unmapped reads: " + gfc.getNumberOfPartiallyUnmappedReads());
		
		assertEquals( 156 / 100.0,fc.getCount(),0.000001);
		
		assertEquals(2.0,gfc.getTotalCount(),0.000001);
		assertEquals(156 / 100.0,gfc.getMappedReadCount(),0.000001);
		assertEquals(44 / 100.0,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(2,gfc.getNumberOfPartiallyUnmappedReads());
	}
	
	
	@Test
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
		
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(11869);
		record2.setCigarString("39M17400N39M");
		
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		assertEquals(0.0,fc.getCount(),0.00001);
		
		assertEquals(1.0,gfc.getTotalCount(),0.000001);
		assertEquals(0.0,gfc.getMappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(1,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(0,gfc.getNumberOfPartiallyUnmappedReads());
	}
	
	
	
	@Test
	public void testCountOneMapped1() throws IOException {
		
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
		
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(12228);
		record2.setCigarString("39M");
		
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		assertEquals(0.5,fc.getCount(),0.00001);
		
		assertEquals(1.0,gfc.getTotalCount(),0.000001);
		assertEquals(0.5,gfc.getMappedReadCount(),0.000001);
		assertEquals(0.5,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(0,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(1,gfc.getNumberOfPartiallyUnmappedReads());
	}
	
	@Test
	public void testCountOneMapped2() throws IOException {
		
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record = new SAMRecord(new SAMFileHeader());
		
		record.setReferenceName("chr1");
		record.setAlignmentStart(12228);
		record.setCigarString("39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record);
		
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		
		record2.setReferenceName("chr1");
		record2.setAlignmentStart( 11869);
		record2.setCigarString("39M");
		
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		FeatureCount fc = gfc.getCounts("ENSG00000223972_chr1_11869_12227");
		assertEquals(0.5,fc.getCount(),0.00001);
		
		assertEquals(1.0,gfc.getTotalCount(),0.000001);
		assertEquals(0.5,gfc.getMappedReadCount(),0.000001);
		assertEquals(0.5,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(0,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(1,gfc.getNumberOfPartiallyUnmappedReads());
	}
	
	@Test
	public void testCountPairFiltered() throws IOException {
		logger.info("testCountPair1");
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		//File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.addFilter(new MaxAlignmentDistanceFilter(5));
		
		//gfc.buildFeatures(gtf, "exon");
		
		
		SAMRecord record1 = new SAMRecord(new SAMFileHeader());
		record1.setReferenceName("chr1");
		record1.setAlignmentStart(11869);
		record1.setCigarString("39M");
		record1.setAttribute("nM", 6);
		
		SAMRecord record2 = new SAMRecord(new SAMFileHeader());
		record2.setReferenceName("chr1");
		record2.setAlignmentStart(12000);
		record2.setCigarString("39M");
		
		SAMRecordPair pair = new SAMRecordPair();
		pair.setMate1(record1);
		pair.setMate2(record2);
		
		gfc.count(pair);
		
		assertEquals(1.0,gfc.getTotalCount(),0.000001);
		assertEquals(0.0,gfc.getMappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getUnmappedReadCount(),0.000001);
		assertEquals(0.0,gfc.getAmbiguousReadCount(),0.000001);
		assertEquals(0,gfc.getNumberOfPartiallyUnmappedReads());
		assertEquals(1,gfc.getNumberOfFilteredReads());
	}
	
	
	@Test
	public void computeGeneLevelExpressionTest() {
		//static List<Feature> computeGeneLevelExpression(Map<String,FeatureCount> featureCounts, Map<String,Feature> geneInfos, double numberOfReads)
		
		Feature f1 = new Feature("chr1","havana","exon",Location.fromBio(100, 200, '+'),0.0,0,"gene_id \"gene_1\";");
		Feature f2 = new Feature("chr1","havana","exon",Location.fromBio(300, 400, '+'),0.0,0,"gene_id \"gene_1\";");
		Feature f3 = new Feature("chr1","havana","exon",Location.fromBio(500, 600, '+'),0.0,0,"gene_id \"gene_1\";");
		Feature f4 = new Feature("chr1","havana","exon",Location.fromBio(700, 800, '+'),0.0,0,"gene_id \"gene_1\";");
		
		FeatureCount fc1 = new FeatureCount(f1);
		fc1.addToCount(100);
		FeatureCount fc2 = new FeatureCount(f2);
		fc2.addToCount(200);
		FeatureCount fc3 = new FeatureCount(f3);
		fc3.addToCount(400);
		FeatureCount fc4 = new FeatureCount(f4);
		fc4.addToCount(500);
		
		Feature g1 = new Feature("chr1","havana","gene",Location.fromBio(100, 1000, '+'),0.0,0,"gene_id \"gene_1\";");
		
		HashMap<String,FeatureCount> fcs = new HashMap<>();
		fcs.put(fc1.getId(), fc1);
		fcs.put(fc2.getId(), fc2);
		fcs.put(fc3.getId(), fc3);
		fcs.put(fc4.getId(), fc4);
		
		HashMap<String, Feature> geneInfos = new HashMap<>();
		geneInfos.put(g1.getAttribute("gene_id"), g1);
		
		double numberOfReads = 1000000;
		
		List<Feature> geneCounts = GTExFeatureCounter.computeGeneLevelExpression(fcs, geneInfos, numberOfReads);
		
		assertEquals(1,geneCounts.size());
		
		Feature gene = geneCounts.get(0);
		logger.info(GTFFeatureRenderer.render(gene));
		
		assertEquals("404", gene.getAttribute("length"));
		assertEquals("1200.0", gene.getAttribute("reads"));
		assertEquals(Double.toString( 1200.0 / 404.0 / numberOfReads * Math.pow(10, 9) ), gene.getAttribute("RPKM"));
		assertEquals("gene_1_chr1_100_200,gene_1_chr1_300_400,gene_1_chr1_500_600,gene_1_chr1_700_800", gene.getAttribute("exons"));

	}
	
	@Test
	public void testStrandOnCounts() throws IOException {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		File gtf = r.getFile();
		
		FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		GTExFeatureCounter gfc = new GTExFeatureCounter(findOverlappingFeatures, new LoggerService());
		gfc.buildFeatures(gtf, "exon");
		
		List<Feature> features = gfc.getCounts();
		
		int found = 0;
		for(Feature feature : features) {
			String geneId = feature.getAttribute("gene_id");
			
			assertTrue(feature.hasAttribute("tss"));
			
			if(geneId.equals("ENSG00000187583")) {
				assertEquals(901877, Integer.parseInt(feature.getAttribute("tss")));
				found++;
			}
			
			
			if(geneId.equals("ENSG00000187634")) {
				assertEquals(860260, Integer.parseInt(feature.getAttribute("tss")));
				found++;
			}
			
			if(geneId.equals("ENSG00000268179")) {
				assertEquals(866445, Integer.parseInt(feature.getAttribute("tss")));
				found++;
			}
			
			if(geneId.equals("ENSG00000236601")) {
				assertEquals(460480, Integer.parseInt(feature.getAttribute("tss")));
				found++;
			}
		}
		
		assertTrue(found >= 4);
		
		
	}
}
