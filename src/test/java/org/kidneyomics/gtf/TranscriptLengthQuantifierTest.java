package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class TranscriptLengthQuantifierTest {

	Logger logger = LoggerFactory.getLogger(TranscriptLengthQuantifierTest.class);
	
	@Test
	public void test() throws IOException {
		TranscriptLengthQuantifier t = TranscriptLengthQuantifier.getTranscriptLengthQuantifier("exon");
		
		Resource r = new ClassPathResource("gencode.head.gz");
		
		GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(r.getFile());
		
		List<Feature> list = reader.readAllLines();
		
		List<Feature> transcripts = new LinkedList<Feature>();
		
		for(Feature f : list) {
			if(f.type().equals("exon")) {
				t.addExon(f);
			} else if(f.type().equals("transcript")) {
				transcripts.add(f);
			}
		}
		
		for(Feature f : transcripts) {
			String id = f.getAttribute("transcript_id");
			logger.info(id + ": " + t.length(id));
		}
		
		// gzcat gencode.head.gz | grep ENST00000609139
		// chr1	HAVANA	exon	763049	763155
		// chr1	HAVANA	exon	764383	764484
		// chr1	HAVANA	exon	776580	776753
		// chr1	HAVANA	exon	787307	787378
		
		int actualLength = (763155 - 763049 + 1) + (764484 - 764383 + 1) + (776753 - 776580 + 1) + (787378 - 787307 + 1);
		
		int calcLength = t.length("ENST00000609139");
		
		assertEquals(actualLength,calcLength);
	}

}
