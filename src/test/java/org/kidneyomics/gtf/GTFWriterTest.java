package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class GTFWriterTest {

	Logger logger = LoggerFactory.getLogger(GTFWriterTest.class);
	
	@Test
	public void test() throws IOException {
		Resource r = new ClassPathResource("gencode.head.gz");
		
		GTFReader reader = GTFReader.getGTFByFile(r.getFile());
		
		List<Feature> list = reader.readAllLines();
		
		File tmpDir = FileUtils.getTempDirectory();
		File tmp = new File(tmpDir.getAbsolutePath() + "/test.gtf");
		
		logger.info("Writing to " + tmp.getAbsolutePath());
		GTFWriter writer = GTFWriter.getGTFWriterForFile(tmp);
		
		writer.write(list);
		writer.close();
		
		
		
		GTFReader reader2 = GTFReader.getGTFByFile(tmp);
		
		List<Feature> list2 = reader2.readAllLines();
		
		for(Feature feature : list2) {
			logger.info(GTFFeatureRenderer.render(feature));
		}
		
		assertEquals(list.size(),list2.size());
		
		if(tmp.exists()) {
			tmp.delete();
		}
		
		if(tmpDir.exists()) {
			tmpDir.delete();
		}
		
	}

}
