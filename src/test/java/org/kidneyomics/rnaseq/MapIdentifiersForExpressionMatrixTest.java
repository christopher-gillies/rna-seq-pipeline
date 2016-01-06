package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class MapIdentifiersForExpressionMatrixTest {

	@Test
	public void test() throws Exception {
		LoggerService service = new LoggerService();
		Logger logger = service.getLogger(this);
		
		ApplicationOptions options = new ApplicationOptions(service);
		
		Resource r = new ClassPathResource("test_genes.txt");
		
		List<String> lines = new LinkedList<>();
		lines.add("9\tAAAAA");
		lines.add("10\tBBBB");
		lines.add("134\tCCCC");
		
		String tmpDir = FileUtils.getTempDirectoryPath();
		File in = new File(tmpDir + "/map.txt");
		FileUtils.writeLines(in, lines);
		
		File out = new File(tmpDir + "/outMatrix.txt");
		
		
		options.setExpressionMatrix(r.getFile().getAbsolutePath());
		options.setFileIn(in.getAbsolutePath());
		options.setFileOut(out.getAbsolutePath());
		
		MapIdentifiersForExpressionMatrix mapper = new MapIdentifiersForExpressionMatrix(service, options);
		mapper.doWork();
		
		
		assertTrue(out.exists());
		
		List<String> outLines = FileUtils.readLines(out);
		assertEquals(10,outLines.size());
		String newHeader = outLines.get(0);
		
		//counting starts from 0
		String[] newHeaderCols = newHeader.split("\t");
		assertEquals(135,newHeaderCols.length);
		
		assertEquals("AAAAA", newHeaderCols[9]);
		assertEquals("BBBB", newHeaderCols[10]);
		assertEquals("CCCC", newHeaderCols[134]);
		
		for(String line : outLines) {
			System.err.println(line);
			String[] cols = line.split("\t");
			assertEquals(135,cols.length);
		}
		
		//Clean up
		in.delete();
		
	}

}
