package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.springframework.stereotype.Component;

@Component
public class SampleGTFReader implements SupportFileReader<SampleGTF> {

	@Override
	public List<SampleGTF> readFile(File f) throws IOException {
		/*
		 * Read sample information
		 */
		List<String> sampleLines = FileUtils.readLines(f);
		List<SampleGTF> sampleGtfs = new ArrayList<SampleGTF>(sampleLines.size());
		for(String line : sampleLines) {
			SampleGTF sample = new SampleGTF(line);
			sampleGtfs.add(sample);
		}
		
		return sampleGtfs;
	}

	@Override
	public List<String> getIds(List<SampleGTF> items) {
		List<String> sampleIds = new ArrayList<String>(items.size());
		for(SampleGTF sg : items) {
			sampleIds.add(sg.getSampleId());
		}
		return sampleIds;
	}

}
