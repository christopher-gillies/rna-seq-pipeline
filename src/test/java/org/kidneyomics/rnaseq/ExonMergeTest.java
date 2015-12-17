package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GTFWriter;
import org.kidneyomics.gtf.GeneOrExonFilter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class ExonMergeTest {

	
	Logger logger = LoggerFactory.getLogger(ExonMergeTest.class);
	LoggerService loggerService = new LoggerService();
	
	public List<Feature> setup(Resource r, File tmpDir, File tmpSample1, File tmpSample2, File infile) throws IOException {
		
		GTFReader reader = GTFReader.getGTFByFile(r.getFile());
		reader.addFilter(new GeneOrExonFilter());
		
		List<Feature> exonsAndGenes = reader.readAllLines();
		
		
		/*
		 * Create sample 1
		 */
		List<Feature> sample1 = new LinkedList<Feature>();
		
		int evenCount = 0;
		for(Feature feature : exonsAndGenes) {
			
			Map<String,String> counts = new HashMap<String,String>();
			
			counts.put("reads", Integer.toString(evenCount));
			counts.put("RPKM", Integer.toString(100000 + evenCount));
			counts.put("tss", Integer.toString(feature.location().bioStart()));
			
			if(feature.type().equals("exon")) {
				counts.put("id", FeatureCount.featureIdMaker(feature) + "_" + feature.getAttribute("transcript_id") );
			} else if(feature.type().equals("gene")) {
				counts.put("length", Integer.toString(feature.location().length()));
			} else {
				throw new IllegalStateException("");
			}
			
			Feature newFeature = GTFFeatureBuilder.addAttributesToFeature(feature, counts);
			
			sample1.add(newFeature);
			evenCount += 2;
		}
		
		
		/*
		 * Create sample 2
		 */
		List<Feature> sample2 = new LinkedList<Feature>();
		
		int oddCount = 1;
		for(Feature feature : exonsAndGenes) {
			
			Map<String,String> counts = new HashMap<String,String>();
			
			counts.put("reads", Integer.toString(oddCount));
			counts.put("RPKM", Integer.toString(100000 + oddCount));
			counts.put("tss", Integer.toString(feature.location().bioStart()));
			
			if(feature.type().equals("exon")) {
				counts.put("id", FeatureCount.featureIdMaker(feature) + "_" + feature.getAttribute("transcript_id"));
			} else if(feature.type().equals("gene")) {
				counts.put("length", Integer.toString(feature.location().length()));
			} else {
				throw new IllegalStateException("");
			}
			
			Feature newFeature = GTFFeatureBuilder.addAttributesToFeature(feature, counts);
			
			sample2.add(newFeature);
			oddCount += 2;
		}
		
		

		
		GTFWriter writer1 = GTFWriter.getGTFWriterForFile(tmpSample1);
		writer1.write(sample1);
		writer1.close();
		
		//logger.info("Sample 1");
		//GTFReader readerSample1 = GTFReader.getGTFByFile(tmpSample1);
		//for(Feature feature : readerSample1) {
		//	logger.info( GTFFeatureRenderer.render(feature));
		//}
		
		
		GTFWriter writer2 = GTFWriter.getGTFWriterForFile(tmpSample2);
		writer2.write(sample2);
		writer2.close();
		
		//logger.info("Sample 2");
		//GTFReader readerSample2 = GTFReader.getGTFByFile(tmpSample2);
		//for(Feature feature : readerSample2) {
		//	logger.info( GTFFeatureRenderer.render(feature));
		//}
		
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
		
		return exonsAndGenes;
		
	}
	
	@Test
	public void test() throws IOException {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmpSample1 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample1.gtf");
		File tmpSample2 = new File(tmpDir.getAbsolutePath() + "/test.tmpSample2.gtf");
		File infile = new File(tmpDir.getAbsolutePath() + "/" + "test.samplelist.txt");
		String outdir = tmpDir.getAbsolutePath();
		File annotationFile = new File(tmpDir.getAbsolutePath() + "/" + "anno.gtf.gz");
		FileUtils.copyFile(r.getFile(), annotationFile);
		
		List<Feature> genesAndExons = setup(r, tmpDir, tmpSample1, tmpSample2, infile);
		
		int exonCount = 0;
		int geneCount = 0;
		for(Feature f : genesAndExons) {
			if(f.type().equals("gene")) {
				geneCount++;
			}
			
			if(f.type().equals("exon")) {
				exonCount++;
			}
		}
		
		ApplicationOptions ao = new ApplicationOptions(new LoggerService());
		ao.setFileIn(infile.getAbsolutePath());
		ao.setOutputDirectory(outdir);
		
		ExonMerge exm = new ExonMerge(new LoggerService(), ao, new SampleGTFReader(), new QuantificationFactory());
		
		exm.writeOutMatrices();
		
		
		for(Map.Entry<String, MutableQuantification> entry : exm.exonCounts.entrySet()) {
			MutableQuantification mq = entry.getValue();
			assertTrue(mq.getSampleExpression("sample1") % 2 == 0);
			assertTrue(mq.getSampleExpression("sample2") % 2 == 1);
		}
		
		for(Map.Entry<String, MutableQuantification> entry : exm.exonRpkm.entrySet()) {
			MutableQuantification mq = entry.getValue();
			assertTrue(mq.getSampleExpression("sample1") % 2 == 0);
			assertTrue(mq.getSampleExpression("sample2") % 2 == 1);
			assertTrue(mq.getSampleExpression("sample1") >= 100000);
			assertTrue(mq.getSampleExpression("sample2") >= 100000);
		}
		
		for(Map.Entry<String, MutableQuantification> entry : exm.geneCounts.entrySet()) {
			MutableQuantification mq = entry.getValue();
			assertTrue(mq.getSampleExpression("sample1") % 2 == 0);
			assertTrue(mq.getSampleExpression("sample2") % 2 == 1);
		}
		
		for(Map.Entry<String, MutableQuantification> entry : exm.geneRpkm.entrySet()) {
			MutableQuantification mq = entry.getValue();
			assertTrue(mq.getSampleExpression("sample1") % 2 == 0);
			assertTrue(mq.getSampleExpression("sample2") % 2 == 1);
			assertTrue(mq.getSampleExpression("sample1") >= 100000);
			assertTrue(mq.getSampleExpression("sample2") >= 100000);
		}
		
		/*
		 * 		
		//write results
		writeResults(exonCounts,outDir + "/exon.counts.txt");
		writeResults(exonRpkm,outDir + "/exon.rpkm.txt");
		writeResults(geneCounts,outDir + "/gene.counts.txt");
		writeResults(geneRpkm,outDir + "/gene.rpkm.txt");
		 */
		File exonCounts = new File(outdir + "/exon.counts.txt");
		File exonRpkm = new File(outdir + "/exon.rpkm.txt");
		File geneCounts = new File(outdir + "/gene.counts.txt");
		File geneRpkm = new File(outdir + "/gene.rpkm.txt");
		
		assertTrue(exonCounts.exists());
		assertEquals(exonCount + 1,FileUtils.readLines(exonCounts).size());
		assertTrue(exonRpkm.exists());
		assertEquals(exonCount + 1,FileUtils.readLines(exonRpkm).size());
		assertTrue(geneCounts.exists());
		assertEquals(geneCount + 1,FileUtils.readLines(geneCounts).size());
		assertTrue(geneRpkm.exists());
		assertEquals(geneCount + 1,FileUtils.readLines(geneRpkm).size());
		
		logger.info("wrote file: " + exonCounts.getAbsolutePath());
		logger.info("wrote file: " + exonRpkm.getAbsolutePath());
		logger.info("wrote file: " + geneCounts.getAbsolutePath());
		logger.info("wrote file: " + geneRpkm.getAbsolutePath());
		
		annotationFile.delete();
		infile.delete();
		tmpSample1.delete();
		tmpSample2.delete();
		exonCounts.delete();
		exonRpkm.delete();
		geneCounts.delete();
		geneRpkm.delete();
		tmpDir.delete();
		
	}

}
