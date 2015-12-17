package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ExonOrGeneQuantificationTest {

	
	Logger logger = LoggerFactory.getLogger(ExonOrGeneQuantificationTest.class);
	
	@Test
	public void testHeaderExon() {
		
		String f = "chr1	ENSEMBL	exon	818043	819983	.	+	.	gene_id \"ENSG00000269308.1\"; id \"ENSG00000269308.1_chr1_818043_819983\"; tss \"818043\"; transcript_id \"ENSG00000269308.1\"; gene_type \"protein_coding\"; gene_status \"NOVEL\"; gene_name \"AL645608.2\"; transcript_type \"protein_coding\"; transcript_status \"NOVEL\"; transcript_name \"AL645608.2\"; level 3;";
		Feature feature = GTFFeatureBuilder.createFromLine(f);
		
		List<String> ids = new LinkedList<>();
		ids.add("sample_1");
		ids.add("sample_2");
		
		ExonOrGeneQuantificationResult result =  ExonOrGeneQuantificationResult.getExonCountsQuantification(feature, ids);
		
		String exp = "exon_id\tgene_id\tgene_name\tgene_type\tchr\ttranscription_start_site\tstart\tend\tlength\tstrand\tsample_1\tsample_2";
		assertEquals(exp,result.printHeader());
		
	}
	
	@Test
	public void testToStringExon() {
		
		String f = "chr1	ENSEMBL	exon	818043	819983	.	+	.	gene_id \"ENSG00000269308.1\"; id \"ENSG00000269308.1_chr1_818043_819983\"; tss \"818043\"; transcript_id \"ENSG00000269308.1\"; gene_type \"protein_coding\"; gene_status \"NOVEL\"; gene_name \"AL645608.2\"; transcript_type \"protein_coding\"; transcript_status \"NOVEL\"; transcript_name \"AL645608.2\"; level 3; reads \"11\"";
		Feature feature = GTFFeatureBuilder.createFromLine(f);
		
		logger.info("feature id: " + feature.getAttribute("id"));
		
		List<String> ids = new LinkedList<>();
		ids.add("sample_1");
		ids.add("sample_2");
		
		ExonOrGeneQuantificationResult result =  ExonOrGeneQuantificationResult.getExonCountsQuantification(feature, ids);
		
		result.putSampleCount("sample_1", feature);
		result.putSampleExpression("sample_2", 13);
		
		String exp = "ENSG00000269308.1_chr1_818043_819983\tENSG00000269308.1\tAL645608.2\tprotein_coding\tchr1\t818043\t818043\t819983\t1941\t+\t11.0\t13.0";
		logger.info(exp);
		logger.info(result.toString());
		assertEquals(exp,result.toString());
		
	}
	
	@Test
	public void testToStringGene() {
		
		String f = "chr1	ENSEMBL	gene	818043	819983	.	+	.	gene_id \"ENSG00000269308.1\"; length \"2000\" tss \"818043\"; transcript_id \"ENSG00000269308.1\"; gene_type \"protein_coding\"; gene_status \"NOVEL\"; gene_name \"AL645608.2\"; transcript_type \"protein_coding\"; transcript_status \"NOVEL\"; transcript_name \"AL645608.2\"; level 3; reads \"11\"";
		Feature feature = GTFFeatureBuilder.createFromLine(f);
		
		//logger.info("feature id: " + feature.getAttribute("id"));
		
		List<String> ids = new LinkedList<>();
		ids.add("sample_1");
		ids.add("sample_2");
		
		ExonOrGeneQuantificationResult result =  ExonOrGeneQuantificationResult.getGeneCountsQuantification(feature, ids);
		
		result.putSampleCount("sample_1", feature);
		result.putSampleExpression("sample_2", 13);
		
		String exp = "ENSG00000269308.1\tAL645608.2\tprotein_coding\tchr1\t818043\t818043\t819983\t2000\t+\t11.0\t13.0";
		logger.info(exp);
		logger.info(result.toString());
		assertEquals(exp,result.toString());
		
	}
	
	
	
	@Test
	public void testHeaderGene() {
		
		String f = "chr1	ENSEMBL	gene	818043	819983	.	+	.	gene_id \"ENSG00000269308.1\"; length \"100\"; tss \"818043\"; transcript_id \"ENSG00000269308.1\"; gene_type \"protein_coding\"; gene_status \"NOVEL\"; gene_name \"AL645608.2\"; transcript_type \"protein_coding\"; transcript_status \"NOVEL\"; transcript_name \"AL645608.2\"; level 3;";
		Feature feature = GTFFeatureBuilder.createFromLine(f);
		
		List<String> ids = new LinkedList<>();
		ids.add("sample_1");
		ids.add("sample_2");
		
		ExonOrGeneQuantificationResult result =  ExonOrGeneQuantificationResult.getGeneCountsQuantification(feature, ids);

		
		String exp = "gene_id\tgene_name\tgene_type\tchr\ttranscription_start_site\tstart\tend\tlength\tstrand\tsample_1\tsample_2";
		assertEquals(exp,result.printHeader());
		
	}

}
