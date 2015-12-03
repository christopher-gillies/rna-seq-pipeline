package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class TranscriptLengthQuantifierTest {

	@Test
	public void test() throws IOException {
		TranscriptLengthQuantifier t = TranscriptLengthQuantifier.getTranscriptLengthQuantifier("exon");
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		GTFReader reader = GTFReader.getGTFByFile(r.getFile());
		
		List<Feature> list = reader.readAllLines();
		
		for(Feature f : list) {
			t.addExon(f);
		}
		
		
		// gzcat gencode.head.gz | grep ENST00000609139

		
		
		
	}

}
