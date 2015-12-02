package org.kidneyomics.gtf;

import static org.junit.Assert.*;


import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GTFFeatureBuilderTest {

	Logger logger = LoggerFactory.getLogger(GTFFeatureBuilderTest.class);
	@Test
	public void test() {
		
		
		String input = "chr1	HAVANA	CDS	906066	906138	.	+	0	gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"PLEKHN1-004\"; exon_number 5;  exon_id \"ENSE00001385460.1\";  level 2; tag \"basic\"; tag \"appris_candidate\"; tag \"CCDS\"; ccdsid \"CCDS53256.1\"; havana_gene \"OTTHUMG00000040756.4\"; havana_transcript \"OTTHUMT00000473255.1\";";
		
		Feature feature = GTFFeatureBuilder.createFromLine(input);
		
		logger.info(feature.getAttribute("tag"));
		
		assertEquals(feature.getAttribute("gene_name"),"PLEKHN1");
		assertEquals(feature.getAttribute("gene_id"),"ENSG00000187583.6");
		assertEquals(feature.getAttribute("gene_type"),"protein_coding");
		assertEquals(feature.getAttribute("transcript_id"),"ENST00000379407.3");
		assertEquals(feature.getAttribute("transcript_type"),"protein_coding");
		
		assertEquals(feature.seqname(),"chr1");
		assertEquals(feature.source(),"HAVANA");
		assertEquals(feature.type(),"CDS");
		assertEquals(feature.frame(),0);
		assertEquals(feature.location().length(),906138 - 906066 + 1);
		assertEquals(feature.location().bioStart(),906066);
		assertEquals(feature.location().bioEnd(),906138);
		assertEquals(feature.location().bioStrand(),'+');
	}

	@Test
	public void test2() {
		
		
		String input = "chr1	HAVANA	CDS	906066	906138	.	-	0	gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"PLEKHN1-004\"; exon_number 5;  exon_id \"ENSE00001385460.1\";  level 2; tag \"basic\"; tag \"appris_candidate\"; tag \"CCDS\"; ccdsid \"CCDS53256.1\"; havana_gene \"OTTHUMG00000040756.4\"; havana_transcript \"OTTHUMT00000473255.1\";";
		
		Feature feature = GTFFeatureBuilder.createFromLine(input);
		
		logger.info(feature.getAttribute("tag"));
		
		assertEquals(feature.getAttribute("gene_name"),"PLEKHN1");
		assertEquals(feature.getAttribute("gene_id"),"ENSG00000187583.6");
		assertEquals(feature.getAttribute("gene_type"),"protein_coding");
		assertEquals(feature.getAttribute("transcript_id"),"ENST00000379407.3");
		assertEquals(feature.getAttribute("transcript_type"),"protein_coding");
		
		assertEquals(feature.seqname(),"chr1");
		assertEquals(feature.source(),"HAVANA");
		assertEquals(feature.type(),"CDS");
		assertEquals(feature.frame(),0);
		assertEquals(feature.location().length(),906138 - 906066 + 1);
		assertEquals(feature.location().bioStart(),906066);
		assertEquals(feature.location().bioEnd(),906138);
		assertEquals(feature.location().bioStrand(),'-');
	}
}
