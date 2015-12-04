package org.kidneyomics.rnaseq;

import static org.junit.Assert.assertEquals;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GeneQuantificationTest {

	
	Logger logger = LoggerFactory.getLogger(GeneQuantificationTest.class);
	@Test
	public void test() {
		
		Feature f1 = new Feature("chr1", "a", "transcript", Location.fromBio(100, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		Feature f2 = new Feature("chr1", "a", "transcript", Location.fromBio(50, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379406.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		Feature f3 = new Feature("chr1", "a", "transcript", Location.fromBio(150, 250, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379405.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		
		
		List<String> samples = new LinkedList<String>();
		
		samples.add("sample1");
		samples.add("sample2");
		
		TranscriptQuantification t1 = new TranscriptQuantification(f1, samples);
		
		t1.putSampleExpression("sample1", 25);
		t1.putSampleExpression("sample2", 50);
		
		
		TranscriptQuantification t2 = new TranscriptQuantification(f2, samples);
		
		t2.putSampleExpression("sample1", 25);
		t2.putSampleExpression("sample2", 50);
		
		
		TranscriptQuantification t3 = new TranscriptQuantification(f3, samples);
		
		t3.putSampleExpression("sample1", 25);
		t3.putSampleExpression("sample2", 50);
		
		List<TranscriptQuantification> tqs = new LinkedList<TranscriptQuantification>();
		tqs.add(t1);
		tqs.add(t2);
		tqs.add(t3);
		
		GeneQuantification gq = new GeneQuantification(tqs);
		
		String headerResult = gq.printHeader();
		String headerExp = "gene_id	gene_name	gene_type	chr	transcription_start_site	start	end	length	strand	sample1	sample2";
		
		logger.info(headerResult);
		logger.info(headerExp);
		
		assertEquals(headerExp,headerResult);
		
		String result = gq.toString();
		String expResult = "ENSG00000187583.6	PLEKHN1	protein_coding	chr1	50	50	250	-1	+	75.0	150.0";
		
		logger.info(result);
		logger.info(expResult);
		
		assertEquals(expResult,result);
		
		
	}
}
