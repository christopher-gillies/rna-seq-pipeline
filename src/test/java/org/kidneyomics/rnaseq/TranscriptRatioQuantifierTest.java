package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TranscriptRatioQuantifierTest {

	Logger logger = LoggerFactory.getLogger(TranscriptRatioQuantifierTest.class);
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
		t3.putSampleExpression("sample2", 75);
		
		List<TranscriptQuantification> tqs = new LinkedList<TranscriptQuantification>();
		tqs.add(t1);
		tqs.add(t2);
		tqs.add(t3);
		
		GeneQuantification gq = new GeneQuantification(tqs);
		List<GeneQuantification> gqs = new LinkedList<GeneQuantification>();
		gqs.add(gq);
		
		TranscriptRatioQuantifier trq = new TranscriptRatioQuantifier(gqs, samples);
		
		List<TranscriptQuantification> tqsResult = trq.updateTranscriptRatios();
		
		assertEquals(3,tqsResult.size());
		
		
		String expResult1 = "ENST00000379407.3	ENSG00000187583.6	PLEKHN1	protein_coding	protein_coding	chr1	100	100	200	-1	+\t" + 25.0 / ( 25.0 + 25.0 + 25.0) + "\t" +  50.0 / (50.0 + 50.0 + 75.0);
		// 25 / ( 25 + 25 + 25)
		// 50 / (50 + 50 + 75)
		String actualResult1 = tqsResult.get(0).toString();
		
		logger.info(expResult1);
		logger.info(actualResult1);

		assertEquals(expResult1, actualResult1);
		
		
		String expResult2 = "ENST00000379406.3	ENSG00000187583.6	PLEKHN1	protein_coding	protein_coding	chr1	50	50	200	-1	+\t" + 25.0 / ( 25.0 + 25.0 + 25.0) + "\t" +  50.0 / (50.0 + 50.0 + 75.0);
		String actualResult2 = tqsResult.get(1).toString();
		
		logger.info(expResult2);
		logger.info(actualResult2);

		assertEquals(expResult2, actualResult2);
		
		
		String expResult3 = "ENST00000379405.3	ENSG00000187583.6	PLEKHN1	protein_coding	protein_coding	chr1	150	150	250	-1	+\t" + 25.0 / ( 25.0 + 25.0 + 25.0) + "\t" +  75.0 / (50.0 + 50.0 + 75.0);
		String actualResult3 = tqsResult.get(2).toString();
		
		logger.info(expResult3);
		logger.info(actualResult3);

		assertEquals(expResult3, actualResult3);
	}

}
