package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.util.HashMap;
import java.util.Map;

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
	
	
	@Test
	public void test3() {
		
		
		String input = "chr1	HAVANA	CDS	906066	906138	.	-	0	gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"PLEKHN1-004\"; exon_number 5;  exon_id \"ENSE00001385460.1\";  level 2; tag \"basic\"; tag \"appris_candidate\"; tag \"CCDS\"; ccdsid \"CCDS53256.1\"; havana_gene \"OTTHUMG00000040756.4\"; havana_transcript \"OTTHUMT00000473255.1\";";
		
		Feature feature = GTFFeatureBuilder.createFromLine(input,true);
		
		logger.info(feature.getAttribute("tag"));
		
		assertEquals(feature.getAttribute("gene_name"),"PLEKHN1");
		assertEquals(feature.getAttribute("gene_id"),"ENSG00000187583");
		assertEquals(feature.getAttribute("gene_type"),"protein_coding");
		assertEquals(feature.getAttribute("transcript_id"),"ENST00000379407");
		assertEquals(feature.getAttribute("transcript_type"),"protein_coding");
		assertEquals(feature.getAttribute("exon_id"),"ENSE00001385460");
		
		assertEquals(feature.seqname(),"chr1");
		assertEquals(feature.source(),"HAVANA");
		assertEquals(feature.type(),"CDS");
		assertEquals(feature.frame(),0);
		assertEquals(feature.location().length(),906138 - 906066 + 1);
		assertEquals(feature.location().bioStart(),906066);
		assertEquals(feature.location().bioEnd(),906138);
		assertEquals(feature.location().bioStrand(),'-');
	}
	
	
	@Test
	public void test4() {
		String input = "chr1	HAVANA	CDS	906066	906138	.	-	0	gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\"; transcript_status \"KNOWN\"; transcript_name \"PLEKHN1-004\"; exon_number 5;  exon_id \"ENSE00001385460.1\";  level 2; tag \"basic\"; tag \"appris_candidate\"; tag \"CCDS\"; ccdsid \"CCDS53256.1\"; havana_gene \"OTTHUMG00000040756.4\"; havana_transcript \"OTTHUMT00000473255.1\";";
		
		Feature feature = GTFFeatureBuilder.createFromLine(input,true);
		Map<String,String> newItems = new HashMap<String,String>();
		
		newItems.put("reads", "1000");
		newItems.put("RPKM", "2000");
		
		Feature feature2 = GTFFeatureBuilder.addAttributesToFeature(feature, newItems);
		
		logger.info(input);
		logger.info(GTFFeatureRenderer.render(feature2));
		
		assertEquals("1000",feature2.getAttribute("reads"));
		assertEquals("2000",feature2.getAttribute("RPKM"));
	}
	
	@Test
	public void test5() {
		String input = "1	lincRNA	transcript	763079	788146	.	+	.	transcript_id \"ENST00000416570_ENST00000448975\"; locus_id \"1:762988-794826W\"; gene_id \"ENSG00000228794\"; reads 0.000000; length 612; RPKM 0.000000";
		Feature feature = GTFFeatureBuilder.createFromLine(input,true);
		String transcriptid = feature.getAttribute("transcript_id");
		assertEquals("ENST00000416570",transcriptid);
				
	}
}
