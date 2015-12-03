package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class GTFReaderTest {

	Logger logger = LoggerFactory.getLogger(GTFReaderTest.class);
	
	@Test
	public void test() throws IOException {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		GTFReader reader = GTFReader.getGTFByFile(r.getFile());
		
		List<Feature> list = reader.readAllLines();
		
		reader.close();
		
		logger.info(list.size() + "");
		
		assertEquals(list.size(),995);
		
		
		Feature first = list.get(0);
		
		assertEquals(first.location().bioStart(),11869);
		assertEquals(first.location().bioEnd(),14412);
		assertEquals(first.location().bioStrand(),'+');
		
		assertEquals(first.getAttribute("gene_name"),"DDX11L1");
		assertEquals(first.getAttribute("transcript_id"),"ENSG00000223972.4");
	}
	
	
	@Test
	public void test2() throws IOException {
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(r.getFile());
		
		List<Feature> list = reader.readAllLines();
		
		reader.close();
		
		logger.info(list.size() + "");
		
		assertEquals(list.size(),995);
		
		
		Feature first = list.get(0);
		
		assertEquals(first.location().bioStart(),11869);
		assertEquals(first.location().bioEnd(),14412);
		assertEquals(first.location().bioStrand(),'+');
		
		assertEquals(first.getAttribute("gene_name"),"DDX11L1");
		assertEquals(first.getAttribute("transcript_id"),"ENSG00000223972");
	}

}
