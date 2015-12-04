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

}
