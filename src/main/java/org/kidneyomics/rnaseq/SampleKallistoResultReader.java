package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.springframework.stereotype.Component;

@Component
public class SampleKallistoResultReader implements SupportFileReader<SampleKallistoResult> {

	@Override
	public List<SampleKallistoResult> readFile(File f) throws IOException {
		
		List<String> sampleLines = FileUtils.readLines(f);
		List<SampleKallistoResult> sampleKallistoResults = new ArrayList<SampleKallistoResult>(sampleLines.size());
		for(String line : sampleLines) {
			SampleKallistoResult sample = new SampleKallistoResult(line);
			sampleKallistoResults.add(sample);
		}
		
		return sampleKallistoResults;
	}

	@Override
	public List<String> getIds(List<SampleKallistoResult> items) {
		List<String> sampleIds = new ArrayList<String>(items.size());
		for(SampleKallistoResult sg : items) {
			sampleIds.add(sg.getSampleId());
		}
		return sampleIds;
	}
	
}
