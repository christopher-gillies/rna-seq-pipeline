package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class KallistoMergeTest {

	@Test
	public void test() throws Exception {
		Resource resource = new ClassPathResource("gencode.head.gz");
		
		//create list file
		//25856	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/KALLISTO//25856//net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/KALLISTO//25856//abundance.txt	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/KALLISTO//25856//abundance.hd5
		//25858	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/KALLISTO//25858//net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/KALLISTO//25858//abundance.txt	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/KALLISTO//25858//abundance.hd5
		
		
		File tmpDir = FileUtils.getTempDirectory();
		
		File listFile = new File(tmpDir.getAbsoluteFile() + "/list.txt");
		File outAbundance1 = new File(tmpDir.getAbsoluteFile() + "/sample1.tsv");
		File outAbundance2 = new File(tmpDir.getAbsoluteFile() + "/sample2.tsv");
		File out = new File(tmpDir.getAbsoluteFile() + "/out.txt");
		
		if(out.exists()) {
			out.delete();
		}
		
		StringBuilder listFileBuilder = new StringBuilder();
		
		listFileBuilder.append("sample1");
		listFileBuilder.append("\t");
		listFileBuilder.append(".");
		listFileBuilder.append("\t");
		listFileBuilder.append(outAbundance1);
		listFileBuilder.append("\t");
		listFileBuilder.append(".");
		listFileBuilder.append("\n");
		
		listFileBuilder.append("sample2");
		listFileBuilder.append("\t");
		listFileBuilder.append(".");
		listFileBuilder.append("\t");
		listFileBuilder.append(outAbundance2);
		listFileBuilder.append("\t");
		listFileBuilder.append(".");
		listFileBuilder.append("\n");
		
		FileUtils.write(listFile, listFileBuilder.toString());
		
		String data1 = "target_id	length	eff_length	est_counts	tpm\n" + 
				"ENST00000456328	8	3	0	1\n" + 
				"ENST00000515242	9	3	0	2\n" + 
				"ENST00000518655	13	5.16667	0	3\n" + 
				"ENST00000450305	23	5.7	0	4\n" + 
				"ENST00000438504	19	5.46667	0	5";
		
		FileUtils.write(outAbundance1, data1);
		
		String data2 = "target_id	length	eff_length	est_counts	tpm\n" + 
				"ENST00000456328	8	4	0	6\n" + 
				"ENST00000515242	9	4	0	7\n" + 
				"ENST00000518655	13	9.16667	0	8\n" + 
				"ENST00000450305	23	9.7	0	9\n" + 
				"ENST00000438504	19	9.46667	0	10";
		
		FileUtils.write(outAbundance2, data2);
		
		/*
		 		String gtfList = applicationOptions.getFileIn();
				String annotation = applicationOptions.getGtf();
				String outmatrix = applicationOptions.getFileOut();
				boolean outCounts = applicationOptions.isOutCounts();
		 */
		LoggerService service = new LoggerService();
		Logger logger = service.getLogger(this);
		
		ApplicationOptions applicationOptions = new ApplicationOptions(service);
		
		applicationOptions.setFileIn(listFile.getAbsolutePath());
		applicationOptions.setGtf(resource.getFile().getAbsolutePath());
		applicationOptions.setFileOut(out.getAbsolutePath());
		applicationOptions.setOutCounts(false);
		
		KallistoMerge merger = new KallistoMerge(service, applicationOptions, new SampleKallistoResultReader());
		
		merger.writeTranscriptMatrix();
		
		assertTrue(out.exists());
		
		List<String> resultLines = FileUtils.readLines(out);
		for(String line : resultLines) {
			logger.info(line);
		}
		
		assertEquals("ENST00000456328	ENSG00000223972	DDX11L1	pseudogene	processed_transcript	chr1	11869	11869	14409	1657	+	1.0	6.0",resultLines.get(1));
		
	}

}
