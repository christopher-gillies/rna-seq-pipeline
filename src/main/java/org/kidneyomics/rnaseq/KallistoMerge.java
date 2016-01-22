package org.kidneyomics.rnaseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.GTFReader;
import org.parboiled.common.StringUtils;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component()
public class KallistoMerge extends TranscriptResultsMerger {

	private final SampleKallistoResultReader reader;
	
	@Autowired
	KallistoMerge(LoggerService loggerService, ApplicationOptions applicationOptions, SampleKallistoResultReader reader) {
		super(loggerService, applicationOptions);
		this.reader = reader;
	}

	private List<SampleKallistoResult> sampleResults;
	
	@Override
	List<String> getReadFileAndGetSampleIds(String fileIn) throws IOException {
		sampleResults = reader.readFile(new File(fileIn));
		Collections.sort(sampleResults);
		return reader.getIds(sampleResults);
	}

	@Override
	void addTranscriptQuantifications(Map<String, TranscriptQuantification> transcripts, boolean outCounts)
			throws Exception {
		
		/*
		 * Read samples
		 */
		
		for(SampleKallistoResult sample : sampleResults) {
			logger.info("Reading kallisto file for sample " + sample.getSampleId());
			File file = sample.getFile();
			try(BufferedReader reader = Files.newBufferedReader(file.toPath(), Charset.defaultCharset())) {
				String header = reader.readLine();
				String line = null;
				while( (line = reader.readLine()) != null) {
					if(StringUtils.isEmpty(line)) {
						continue;
					}
					
					KallistoResult lineResult = new KallistoResult(line);
					String transcriptId = lineResult.getTranscriptId();
					
					if(transcripts.containsKey(transcriptId)) {
						TranscriptQuantification tq = transcripts.get(transcriptId);
						if(outCounts) {
							tq.putSampleExpression(sample.getSampleId(), lineResult.getCounts());
						} else {
							tq.putSampleExpression(sample.getSampleId(), lineResult.getTpm());
						}
					} else {
						throw new Exception("Invalid transcript id " + transcriptId);
					}
				}
			}
		}
		
	}

}
