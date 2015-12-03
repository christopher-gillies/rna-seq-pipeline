package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.util.StringUtils;

public class FluxMergeTest {

	Logger logger = LoggerFactory.getLogger(FluxMergeTest.class);
	@Test
	public void testReads() throws Exception {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
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
			
			counts.put("reads", Integer.toString(evenCount));
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
			
			counts.put("reads", Integer.toString(oddCount));
			Feature newFeature = GTFFeatureBuilder.addAttributesToFeature(feature, counts);
			
			sample2.add(newFeature);
			oddCount += 2;
		}
		
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		
		
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
		
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		StringBuilder sb = new StringBuilder();
		sb.append("sample1\t");
		sb.append(tmpSample1.getAbsolutePath());
		sb.append("\n");
		sb.append("sample2\t");
		sb.append(tmpSample2.getAbsolutePath());
		sb.append("\n");
		FileUtils.write(infile, sb.toString());
		
		LoggerService loggerService = new LoggerService();
		ApplicationOptions appOpts = new ApplicationOptions(loggerService);
		
		File outfile = new File(tmpDir.getAbsolutePath() + "/" + "test.matrix");
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		
		appOpts.setFileIn(infile.getAbsolutePath());
		appOpts.setFileOut(outfile.getAbsolutePath());
		appOpts.setGtf(annotationFile.getAbsolutePath());
		appOpts.setOutCounts(true);
		
		//FluxMerge
		
		FluxMerge fluxMerge = new FluxMerge(loggerService, appOpts, new SampleGTFReader());
		List<TranscriptQuantification> tqs = fluxMerge.getTranscriptQuantifications(appOpts.getFileIn(), appOpts.getGtf(), appOpts.isOutCounts());
		
		assertEquals(transcripts.size(),tqs.size());
		
		for(TranscriptQuantification tq : tqs) {
			double sample1Exp = tq.getSampleExpression("sample1");
			double sample2Exp = tq.getSampleExpression("sample2");
			
			assertTrue(sample1Exp % 2 == 0);
			assertTrue(sample2Exp % 2 == 1);
		}
		
		assertTrue(tqs.get(0).getFeature().location().bioStart() == 11869);
		
		
		fluxMerge.writeTranscriptMatrix(tqs, outfile.getAbsolutePath());
		List<String> outMatrixTextLines = FileUtils.readLines(outfile);
		
		assertTrue(outMatrixTextLines.size() == tqs.size()  + 1);
		
		assertEquals("transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	start	end	length	strand	sample1	sample2",outMatrixTextLines.get(0));
		
		assertEquals("ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	14409	2541	+	0.0	1.0",outMatrixTextLines.get(1));
		
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

}
