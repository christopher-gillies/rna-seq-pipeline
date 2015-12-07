package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GTFWriter;
import org.kidneyomics.rnaseq.FluxMerge.TranscriptRatioResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.util.StringUtils;

public class FluxMergeTest {

	Logger logger = LoggerFactory.getLogger(FluxMergeTest.class);
	LoggerService loggerService = new LoggerService();
	
	public List<Feature> setup(Resource r, File tmpDir, File tmpSample1, File tmpSample2, File infile , String countType) throws IOException {
		
		GTFReader reader = GTFReader.getGTFByFile(r.getFile());
		
		List<Feature> transcripts = new LinkedList<Feature>();
		
		for(Feature feature : reader) {
			if(feature.type().equals("transcript")) {
				transcripts.add(feature);
			}
		}
		
		
		/*
		 * Create sample 1
		 */
		List<Feature> sample1 = new LinkedList<Feature>();
		
		int evenCount = 0;
		for(Feature feature : transcripts) {
			
			Map<String,String> counts = new HashMap<String,String>();
			
			counts.put(countType, Integer.toString(evenCount));
			Feature newFeature = GTFFeatureBuilder.addAttributesToFeature(feature, counts);
			
			sample1.add(newFeature);
			evenCount += 2;
		}
		
		
		/*
		 * Create sample 2
		 */
		List<Feature> sample2 = new LinkedList<Feature>();
		
		int oddCount = 1;
		for(Feature feature : transcripts) {
			
			Map<String,String> counts = new HashMap<String,String>();
			
			counts.put(countType, Integer.toString(oddCount));
			Feature newFeature = GTFFeatureBuilder.addAttributesToFeature(feature, counts);
			
			sample2.add(newFeature);
			oddCount += 2;
		}
		
		

		
		GTFWriter writer1 = GTFWriter.getGTFWriterForFile(tmpSample1);
		writer1.write(sample1);
		writer1.close();
		
		logger.info("Sample 1");
		GTFReader readerSample1 = GTFReader.getGTFByFile(tmpSample1);
		for(Feature feature : readerSample1) {
			logger.info( GTFFeatureRenderer.render(feature));
		}
		
		
		GTFWriter writer2 = GTFWriter.getGTFWriterForFile(tmpSample2);
		writer2.write(sample2);
		writer2.close();
		
		logger.info("Sample 2");
		GTFReader readerSample2 = GTFReader.getGTFByFile(tmpSample2);
		for(Feature feature : readerSample2) {
			logger.info( GTFFeatureRenderer.render(feature));
		}
		
		/*
		 * Create Dependencies 
		 */
		
		StringBuilder sb = new StringBuilder();
		sb.append("sample1\t");
		sb.append(tmpSample1.getAbsolutePath());
		sb.append("\n");
		sb.append("sample2\t");
		sb.append(tmpSample2.getAbsolutePath());
		sb.append("\n");
		FileUtils.write(infile, sb.toString());
		
		return transcripts;
		
	}
	
	@Test
	public void testReads() throws Exception {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		File outfile = new File(tmpDir.getAbsolutePath() + "/" + "test.matrix");
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		List<Feature> transcripts = setup(r, tmpDir, tmpSample1, tmpSample2, infile, "reads");
		
		ApplicationOptions appOpts = new ApplicationOptions(loggerService);
		appOpts.setFileIn(infile.getAbsolutePath());
		appOpts.setFileOut(outfile.getAbsolutePath());
		appOpts.setGtf(annotationFile.getAbsolutePath());
		appOpts.setOutCounts(true);
		
		//FluxMerge
		
		FluxMerge fluxMerge = new FluxMerge(loggerService, appOpts, new SampleGTFReader());
		List<TranscriptQuantification> tqs = fluxMerge.getTranscriptQuantifications(appOpts.getFileIn(), appOpts.getGtf(), appOpts.isOutCounts()).getTranscriptQuantifications();
		
		assertEquals(transcripts.size(),tqs.size());
		
		for(TranscriptQuantification tq : tqs) {
			double sample1Exp = tq.getSampleExpression("sample1");
			double sample2Exp = tq.getSampleExpression("sample2");
			
			assertTrue(sample1Exp % 2 == 0);
			assertTrue(sample2Exp % 2 == 1);
		}
		
		assertTrue(tqs.get(0).getFeature().location().bioStart() == 11869);
		
		
		fluxMerge.writeQuantificationMatrix(tqs, outfile.getAbsolutePath());
		List<String> outMatrixTextLines = FileUtils.readLines(outfile);
		
		assertTrue(outMatrixTextLines.size() == tqs.size()  + 1);
		
		assertEquals("transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2",outMatrixTextLines.get(0));
		
		assertEquals("ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	0.0	1.0",outMatrixTextLines.get(1));
		
		logger.info("\n" + StringUtils.collectionToDelimitedString(outMatrixTextLines, "\n"));
		/*
		 * Delete temporary files
		 */
			
		
		if(infile.exists()) {
			infile.delete();
		}
		
		if(outfile.exists()) {
			outfile.delete();
		}
	
		if(annotationFile.exists()) {
			annotationFile.delete();
		}
		
		if(tmpSample1.exists()) {
			tmpSample1.delete();
		}
		
		if(tmpSample2.exists()) {
			tmpSample2.delete();
		}
		
		if(tmpDir.exists()) {
			tmpDir.delete();
		}
		
		
		
	}
	
	
	@Test
	public void testRPKM() throws Exception {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		File outfile = new File(tmpDir.getAbsolutePath() + "/" + "test.matrix");
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		//set RPKM feature
		List<Feature> transcripts = setup(r, tmpDir, tmpSample1, tmpSample2, infile, "RPKM");
		
		ApplicationOptions appOpts = new ApplicationOptions(loggerService);
		appOpts.setFileIn(infile.getAbsolutePath());
		appOpts.setFileOut(outfile.getAbsolutePath());
		appOpts.setGtf(annotationFile.getAbsolutePath());
		
		//use RPKM
		appOpts.setOutCounts(false);
		
		//FluxMerge
		
		FluxMerge fluxMerge = new FluxMerge(loggerService, appOpts, new SampleGTFReader());
		List<TranscriptQuantification> tqs = fluxMerge.getTranscriptQuantifications(appOpts.getFileIn(), appOpts.getGtf(), appOpts.isOutCounts()).getTranscriptQuantifications();
		
		assertEquals(transcripts.size(),tqs.size());
		
		for(TranscriptQuantification tq : tqs) {
			double sample1Exp = tq.getSampleExpression("sample1");
			double sample2Exp = tq.getSampleExpression("sample2");
			
			assertTrue(sample1Exp % 2 == 0);
			assertTrue(sample2Exp % 2 == 1);
		}
		
		assertTrue(tqs.get(0).getFeature().location().bioStart() == 11869);
		
		
		fluxMerge.writeQuantificationMatrix(tqs, outfile.getAbsolutePath());
		List<String> outMatrixTextLines = FileUtils.readLines(outfile);
		
		assertTrue(outMatrixTextLines.size() == tqs.size()  + 1);
		
		assertEquals("transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2",outMatrixTextLines.get(0));
		
		assertEquals("ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	0.0	1.0",outMatrixTextLines.get(1));
		
		logger.info("\n" + StringUtils.collectionToDelimitedString(outMatrixTextLines, "\n"));
		/*
		 * Delete temporary files
		 */
			
		
		if(infile.exists()) {
			infile.delete();
		}
		
		if(outfile.exists()) {
			outfile.delete();
		}
	
		if(annotationFile.exists()) {
			annotationFile.delete();
		}
		
		if(tmpSample1.exists()) {
			tmpSample1.delete();
		}
		
		if(tmpSample2.exists()) {
			tmpSample2.delete();
		}
		
		if(tmpDir.exists()) {
			tmpDir.delete();
		}
		
		
		
	}
	
	
	@Test
	public void testGeneLevelRPKM() throws Exception {
		//TODO
		/*
		 * 
		 *	transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2
		 *	ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	0.0	1.0
		 *	ENST00000515242	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11872	11872	14412	1653	+	2.0	3.0
		 *	ENST00000518655	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11874	11874	14409	1483	+	4.0	5.0
		 *	ENST00000450305	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	12010	12010	13670	632	+	6.0	7.0
		 *
		 *
		 *	ENST00000423562	ENSG00000227232	WASH7P	pseudogene	unprocessed_pseudogene	chr1	29370	14363	29370	1669	-	12.0	13.0
		 *	ENST00000541675	ENSG00000227232	WASH7P	pseudogene	unprocessed_pseudogene	chr1	24886	14363	24886	1416	-	10.0	11.0
		 *	ENST00000438504	ENSG00000227232	WASH7P	pseudogene	unprocessed_pseudogene	chr1	29370	14363	29370	1783	-	8.0	9.0
		 *	ENST00000488147	ENSG00000227232	WASH7P	pseudogene	unprocessed_pseudogene	chr1	29570	14404	29570	1351	-	14.0	15.0
		 *	ENST00000538476	ENSG00000227232	WASH7P	pseudogene	unprocessed_pseudogene	chr1	29806	14411	29806	1583	-	16.0	17.0
		 */
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		File outfile = new File(tmpDir.getAbsolutePath() + "/" + "test.matrix");
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		//set RPKM feature
		List<Feature> transcripts = setup(r, tmpDir, tmpSample1, tmpSample2, infile, "RPKM");
		
		HashSet<String> geneIds = new HashSet<String>();
		
		for(Feature f : transcripts) {
			String geneId = f.getAttribute("gene_id");
			geneIds.add(geneId);
		}
		
		ApplicationOptions appOpts = new ApplicationOptions(loggerService);
		appOpts.setFileIn(infile.getAbsolutePath());
		appOpts.setFileOut(outfile.getAbsolutePath());
		appOpts.setGtf(annotationFile.getAbsolutePath());
		
		//use RPKM
		appOpts.setOutCounts(false);
		
		//FluxMerge
		
		FluxMerge fluxMerge = new FluxMerge(loggerService, appOpts, new SampleGTFReader());
		fluxMerge.writeGeneMatrix();
		List<String> outMatrixTextLines = FileUtils.readLines(outfile);
		
		assertEquals(geneIds.size() + 1,outMatrixTextLines.size());
		
		logger.info("\n" + StringUtils.collectionToDelimitedString(outMatrixTextLines, "\n"));
		
		assertEquals("gene_id	gene_name	gene_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2",outMatrixTextLines.get(0));
		
		assertEquals("ENSG00000223972	DDX11L1	pseudogene	chr1	11869	11869	14412	1756	+	12.0	16.0",outMatrixTextLines.get(1));
		
		assertEquals("ENSG00000227232	WASH7P	pseudogene	chr1	29806	14363	29806	2073	-	60.0	65.0",outMatrixTextLines.get(2));
		
		
		/*
		 * Delete temporary files
		 */
			
		
		if(infile.exists()) {
			infile.delete();
		}
		
		if(outfile.exists()) {
			outfile.delete();
		}
	
		if(annotationFile.exists()) {
			annotationFile.delete();
		}
		
		if(tmpSample1.exists()) {
			tmpSample1.delete();
		}
		
		if(tmpSample2.exists()) {
			tmpSample2.delete();
		}
		
		if(tmpDir.exists()) {
			tmpDir.delete();
		}
		
	}
	
	
	@Test
	public void testTranscriptRatioRPKM() throws Exception {
		//TODO
		/*
		 * 
		 *	transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2
		 *	ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	0.0	1.0
		 *	ENST00000515242	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11872	11872	14412	1653	+	2.0	3.0
		 *	ENST00000518655	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11874	11874	14409	1483	+	4.0	5.0
		 *	ENST00000450305	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	12010	12010	13670	632	+	6.0	7.0
		 *
		 *
		 */
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		File outfile = new File(tmpDir.getAbsolutePath() + "/" + "test.matrix");
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		//set RPKM feature
		List<Feature> transcripts = setup(r, tmpDir, tmpSample1, tmpSample2, infile, "RPKM");
		
		HashSet<String> geneIds = new HashSet<String>();
		
		for(Feature f : transcripts) {
			String geneId = f.getAttribute("gene_id");
			geneIds.add(geneId);
		}
		
		ApplicationOptions appOpts = new ApplicationOptions(loggerService);
		appOpts.setFileIn(infile.getAbsolutePath());
		appOpts.setFileOut(outfile.getAbsolutePath());
		appOpts.setGtf(annotationFile.getAbsolutePath());
		
		//use RPKM
		appOpts.setOutCounts(false);
		
		//FluxMerge
		
		FluxMerge fluxMerge = new FluxMerge(loggerService, appOpts, new SampleGTFReader());
		fluxMerge.writeTranscriptRatioMatrix();
		List<String> outMatrixTextLines = FileUtils.readLines(outfile);
		
		assertEquals(transcripts.size() + 1,outMatrixTextLines.size());
		
		logger.info("\n" + outMatrixTextLines.get(0));
		logger.info("\n" + outMatrixTextLines.get(1));
		logger.info("\n" + outMatrixTextLines.get(2));
		
		assertEquals("transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2",outMatrixTextLines.get(0));
		
		assertEquals("ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	0.0\t" + (1.0 / (1.0 + 3.0 + 5.0 + 7.0)),outMatrixTextLines.get(1));
		
		assertEquals("ENST00000515242	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11872	11872	14412	1653	+\t" + (2.0 / (2.0 + 4.0 + 6.0)) + "\t" + (3.0 / (1.0 + 3.0 + 5.0 + 7.0)),outMatrixTextLines.get(2));
		
		
		/*
		 * Delete temporary files
		 */
			
		
		if(infile.exists()) {
			infile.delete();
		}
		
		if(outfile.exists()) {
			outfile.delete();
		}
	
		if(annotationFile.exists()) {
			annotationFile.delete();
		}
		
		if(tmpSample1.exists()) {
			tmpSample1.delete();
		}
		
		if(tmpSample2.exists()) {
			tmpSample2.delete();
		}
		
		if(tmpDir.exists()) {
			tmpDir.delete();
		}
		
	}
	
	
	
	@Test
	public void testTranscriptRatioRPKM2() throws Exception {
		//TODO
		/*
		 * 
		 *	transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2
		 *	ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	0.0	1.0
		 *	ENST00000515242	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11872	11872	14412	1653	+	2.0	3.0
		 *	ENST00000518655	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	11874	11874	14409	1483	+	4.0	5.0
		 *	ENST00000450305	ENSG00000223972	DDX11L1	pseudogene	transcribed_unprocessed_pseudogene	chr1	12010	12010	13670	632	+	6.0	7.0
		 *
		 *
		 */
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		//File outfile = new File(tmpDir.getAbsolutePath() + "/" + "test.matrix");
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		//set RPKM feature
		List<Feature> transcripts = setup(r, tmpDir, tmpSample1, tmpSample2, infile, "RPKM");
		
		HashSet<String> geneIds = new HashSet<String>();
		
		for(Feature f : transcripts) {
			String geneId = f.getAttribute("gene_id");
			geneIds.add(geneId);
		}
		
		ApplicationOptions appOpts = new ApplicationOptions(loggerService);
		appOpts.setFileIn(infile.getAbsolutePath());
		//appOpts.setFileOut(outfile.getAbsolutePath());
		appOpts.setGtf(annotationFile.getAbsolutePath());
		
		//use RPKM
		appOpts.setOutCounts(false);
		
		//FluxMerge
		FluxMerge fluxMerge = new FluxMerge(loggerService, appOpts, new SampleGTFReader());
		TranscriptRatioResult trr = fluxMerge.getTranscriptRatios();
		trr.transcriptQuantificationResult.transcriptQuantifications = trr.transcriptRatios;
		List<GeneQuantification> gqs = fluxMerge.getGeneQuantifications(trr.transcriptQuantificationResult);
		
		for(GeneQuantification gq : gqs) {
			for(String id : trr.transcriptQuantificationResult.sampleIds) {
				logger.info("SAMPLE: "+ id + "GENE: " + gq.getGeneId() +  " EXPRESSION: " + gq.getSampleExpression(id));
				assertEquals(1.0,gq.getSampleExpression(id),0.001);
			}
		}
		
		
		/*
		 * Delete temporary files
		 */
			
		
		if(infile.exists()) {
			infile.delete();
		}
		

	
		if(annotationFile.exists()) {
			annotationFile.delete();
		}
		
		if(tmpSample1.exists()) {
			tmpSample1.delete();
		}
		
		if(tmpSample2.exists()) {
			tmpSample2.delete();
		}
		
		if(tmpDir.exists()) {
			tmpDir.delete();
		}
		
	}
	
	
}
