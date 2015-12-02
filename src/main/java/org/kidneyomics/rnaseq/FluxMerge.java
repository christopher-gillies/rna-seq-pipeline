package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.VersionTrimmer;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.apache.commons.io.FileUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;

@Component
public class FluxMerge {
	/**
	 * This class will perform the merging of the individual gtf files into a matrix
	 * Also there will be options for read-level, rpkm, gene-level
	 * This will be built on top of biojava genome parser
	 */
	
	

	ApplicationOptions applicationOptions;
	
	Logger logger;
	
	@Autowired
	FluxMerge(LoggerService loggerService, ApplicationOptions applicationOptions) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
	}
	
	
	public void writeTranscriptMatrix() throws Exception {
		String gtfList = applicationOptions.getFileIn();
		String annotation = applicationOptions.getGtf();
		String outmatrix = applicationOptions.getFileOut();
		boolean outCounts = applicationOptions.isOutCounts();
		
		
		/*
		 * Read sample information
		 */
		List<String> sampleLines = FileUtils.readLines(new File(gtfList));
		List<String> sampleIds = new ArrayList<String>(sampleLines.size());
		List<SampleGTF> sampleGtfs = new ArrayList<SampleGTF>(sampleLines.size());
		for(String line : sampleLines) {
			SampleGTF sample = new SampleGTF(line);
			sampleGtfs.add(sample);
			sampleIds.add(sample.getSampleId());
		}
		
		
		
		/*
		 * Read gencode transcripts
		 */
		HashMap<String,TranscriptQuantification> transcripts = new HashMap<String,TranscriptQuantification>();
		GTFReader reader = GTFReader.getGTFByFile(new File(annotation));
		for(Feature feature : reader) {
			/*
			 * If it is a transcript..
			 */
			if(feature.type().equals("transcript")) {
				String transcriptId = VersionTrimmer.trim(feature.getAttribute("transcript_id"));
				
				transcripts.put(transcriptId, new TranscriptQuantification(feature,sampleIds));
			}
		}
		reader.close();
		
		/*
		 * Read samples
		 */
		
		for(SampleGTF sample : sampleGtfs) {
			GTFReader readerForSample = GTFReader.getGTFByFile(sample.getFile());
			for(Feature feature : readerForSample) {
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
			readerForSample.close();
		}
		
		/*
		 * write results
		 */
		StringBuilder sb = new StringBuilder();
		Collection<TranscriptQuantification> listOfTranscriptQuantifications = transcripts.values();
		
		/*
		 * Perform sort ...
		 */
		
		for(TranscriptQuantification tq : listOfTranscriptQuantifications) {
			tq.toString(sb);
		}
		
		FileUtils.write(new File(outmatrix), sb.toString());
		
		
		
	}
}
