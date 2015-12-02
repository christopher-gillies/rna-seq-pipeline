package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TranscriptQuantificationTest {

	Logger logger = LoggerFactory.getLogger(TranscriptQuantificationTest.class);
	@Test
	public void test() {
		
		Feature f1 = new Feature("chr1", "a", "transcript", Location.fromBio(100, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		
		List<String> samples = new LinkedList<String>();
		
		samples.add("sample1");
		samples.add("sample2");
		
		TranscriptQuantification t = new TranscriptQuantification(f1, samples);
		
		t.putSampleExpression("sample1", 25);
		t.putSampleExpression("sample2", 50);
		
		String headerResult = t.printHeader();
		String headerExp = "transcript_id	gene_id	gene_name	gene_type	transcript_type	chr	start	end	length	strand	sample1	sample2";
		assertEquals(headerExp,headerResult);
		
		String result = t.toString();
		String expResult = "ENST00000379407.3	ENSG00000187583.6	PLEKHN1	protein_coding	protein_coding	chr1	100	200	101	+	25.0	50.0";
		assertEquals(expResult,result);
		
		logger.info(headerResult);
		logger.info(result);
		
	}

}
