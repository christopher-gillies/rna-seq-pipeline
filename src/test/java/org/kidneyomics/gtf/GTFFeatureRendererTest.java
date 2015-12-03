package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GTFFeatureRendererTest {

	Logger logger = LoggerFactory.getLogger(GTFFeatureRendererTest.class);
	
	@Test
	public void test() {
		
		String input = "chr1	HAVANA	CDS	906066	906138	.	+	0	gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"PLEKHN1-004\"; exon_number 5;  exon_id \"ENSE00001385460.1\";  level 2; tag \"basic\"; tag \"appris_candidate\"; tag \"CCDS\"; ccdsid \"CCDS53256.1\"; havana_gene \"OTTHUMG00000040756.4\"; havana_transcript \"OTTHUMT00000473255.1\";";
		
		Feature feature = GTFFeatureBuilder.createFromLine(input);
		
		String results = GTFFeatureRenderer.render(feature);
		
		logger.info(input);
		logger.info(results);
		
		assertEquals(input,results);
	}

}
