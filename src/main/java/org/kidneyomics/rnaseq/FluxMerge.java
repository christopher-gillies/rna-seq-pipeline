package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.GTFReader;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class FluxMerge extends TranscriptResultsMerger {
	/**
	 * This class will perform the merging of the individual gtf files into a matrix
	 * Also there will be options for read-level, rpkm, gene-level
	 * This will be built on top of biojava genome parser
	 */
	
	SampleGTFReader sampleGtfReader;
	
	@Autowired
	FluxMerge(LoggerService loggerService, ApplicationOptions applicationOptions, SampleGTFReader sampleGtfReader) {
		super(loggerService,applicationOptions);
		this.sampleGtfReader = sampleGtfReader;
	}
	
		

	private List<SampleGTF> sampleGtfs;
	
	@Override
	List<String> getReadFileAndGetSampleIds(String fileIn) throws IOException {
		sampleGtfs = sampleGtfReader.readFile(new File(fileIn));
		Collections.sort(sampleGtfs);
		List<String> sampleIds = sampleGtfReader.getIds(sampleGtfs);
		return sampleIds;
	}


	@Override
	void addTranscriptQuantifications(Map<String,TranscriptQuantification> transcripts, boolean outCounts) throws Exception {
		/*
		 * Read samples
		 */
		
		for(SampleGTF sample : sampleGtfs) {
			logger.info("Reading gtf file for sample " + sample.getSampleId());
			GTFReader readerForSample = GTFReader.getGTFByFileNoEnsemblVersion(sample.getFile());
			for(Feature feature : readerForSample) {
				if(feature.type().equals("transcript")) {
					String transcriptId = feature.getAttribute("transcript_id");
					if(transcripts.containsKey(transcriptId)) {
						TranscriptQuantification tq = transcripts.get(transcriptId);
						if(outCounts) {
							tq.putSampleCount(sample.getSampleId(), feature);
						} else {
							tq.putSampleRPKM(sample.getSampleId(), feature);
						}
					} else {
						throw new Exception("Invalid transcript id " + transcriptId);
					}
				}
			}
			readerForSample.close();
		}
	}

	
}


