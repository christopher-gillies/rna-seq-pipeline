package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class GeneLengthQuantifierTest {

	Logger logger = LoggerFactory.getLogger(GeneLengthQuantifierTest.class);
	
	@Test
	public void test() {


		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379406.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(150, 250, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379405.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		

		GeneLengthQuantifier geneLengthQuantifier = GeneLengthQuantifier.getGeneLengthQuantifier("exon");
		
		
		geneLengthQuantifier.addExon(f1);
		geneLengthQuantifier.addExon(f2);
		geneLengthQuantifier.addExon(f3);
		
		int length = geneLengthQuantifier.length("ENSG00000187583.6");
		
		int lengthExp = (250 - 50 + 1);
		
		assertEquals(lengthExp,length);

	}

	@Test
	public void test2() {


		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379407.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 200, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379406.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(230, 250, '+'), 0.0, 0, "gene_id \"ENSG00000187583.6\"; transcript_id \"ENST00000379405.3\"; gene_type \"protein_coding\"; gene_status \"KNOWN\"; gene_name \"PLEKHN1\"; transcript_type \"protein_coding\";");
		

		GeneLengthQuantifier geneLengthQuantifier = GeneLengthQuantifier.getGeneLengthQuantifier("exon");
		
		
		geneLengthQuantifier.addExon(f1);
		geneLengthQuantifier.addExon(f2);
		geneLengthQuantifier.addExon(f3);
		
		int length = geneLengthQuantifier.length("ENSG00000187583.6");
		
		int lengthExp = (200 - 50 + 1) + (250 - 230 + 1);
		
		assertEquals(lengthExp,length);

	}
	
	
	@Test
	public void test3() {

		/*
		 * 
11869	12227
11872	12227
11874	12227

12010	12057

12179	12227

12595	12721

12613	12697

12613	12721

12975	13052
13221	13374
13221	14409

13225	14412
13403	13655
13453	13670
13661	14409
		 */

		Feature f1 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(11869, 12227, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f2 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(11872, 12227, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f3 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(11874, 12227, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f4 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(12010, 12057, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f5 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(12179, 12227, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		
		Feature f6 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(12595, 12721, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f7 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(12613, 12697, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f8 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(12613, 12721, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		
		Feature f9 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(12975, 13052, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		
		Feature f10 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(13221, 13374, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f11 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(13221, 14409, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f12 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(13225, 14412, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f13 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(13403, 13655, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f14 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(13453, 13670, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");
		Feature f15 = new Feature("chr1", "HAVANA", "exon", Location.fromBio(13661, 14409, '+'), 0.0, 0, "gene_id \"ENSG00000223972.4\";");

		
		GeneLengthQuantifier geneLengthQuantifier = GeneLengthQuantifier.getGeneLengthQuantifier("exon");
		
		
		geneLengthQuantifier.addExon(f1);
		geneLengthQuantifier.addExon(f2);
		geneLengthQuantifier.addExon(f3);
		geneLengthQuantifier.addExon(f4);
		geneLengthQuantifier.addExon(f5);
		
		geneLengthQuantifier.addExon(f6);
		geneLengthQuantifier.addExon(f7);
		geneLengthQuantifier.addExon(f8);
		
		geneLengthQuantifier.addExon(f9);
		
		geneLengthQuantifier.addExon(f10);
		geneLengthQuantifier.addExon(f11);
		geneLengthQuantifier.addExon(f12);
		geneLengthQuantifier.addExon(f13);
		geneLengthQuantifier.addExon(f14);
		geneLengthQuantifier.addExon(f15);
		
		int length = geneLengthQuantifier.length("ENSG00000223972.4");
		
		int lengthExp = (12227 - 11869 + 1) + (12721 - 12595 + 1) + (13052 - 12975 + 1) + (14412 - 13221 + 1);
		
		
		assertEquals(lengthExp,length);

	}

}
