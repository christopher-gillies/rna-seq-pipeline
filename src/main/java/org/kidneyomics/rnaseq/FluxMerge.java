package org.kidneyomics.rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
	
	SampleGTFReader sampleGtfReader;
	
	@Autowired
	FluxMerge(LoggerService loggerService, ApplicationOptions applicationOptions, SampleGTFReader sampleGtfReader) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
		this.sampleGtfReader = sampleGtfReader;
	}
	
	
	public void writeTranscriptMatrix() throws Exception {
		String gtfList = applicationOptions.getFileIn();
		String annotation = applicationOptions.getGtf();
		String outmatrix = applicationOptions.getFileOut();
		boolean outCounts = applicationOptions.isOutCounts();
		
		
		Collection<TranscriptQuantification> listOfTranscriptQuantifications = getTranscriptQuantifications(gtfList, annotation, outCounts);
		
		writeTranscriptMatrix(listOfTranscriptQuantifications,outmatrix);
		
	}
	
	protected List<TranscriptQuantification> getTranscriptQuantifications(String gtfList, String annotation, boolean outCounts) throws Exception {
		
		/*
		 * Read sample information
		 */
		List<SampleGTF> sampleGtfs = sampleGtfReader.readFile(new File(gtfList));
		List<String> sampleIds = sampleGtfReader.getIds(sampleGtfs);
		
		
		
		/*
		 * Read gencode transcripts
		 */
		HashMap<String,TranscriptQuantification> transcripts = new HashMap<String,TranscriptQuantification>();
		GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(new File(annotation));
		for(Feature feature : reader) {
			/*
			 * If it is a transcript..
			 */
			if(feature.type().equals("transcript")) {
				String transcriptId = feature.getAttribute("transcript_id");
				
				transcripts.put(transcriptId, new TranscriptQuantification(feature,sampleIds));
			}
		}
		reader.close();
		
		/*
		 * Read samples
		 */
		
		for(SampleGTF sample : sampleGtfs) {
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
		
		/*
		 * return results
		 */
		Collection<TranscriptQuantification> listOfTranscriptQuantifications = transcripts.values();
		List<TranscriptQuantification> list = new ArrayList<TranscriptQuantification>(listOfTranscriptQuantifications);
		logger.info("Sorting results...");
		Collections.sort(list);
		return list;
	}
	
	protected void writeTranscriptMatrix(Collection<TranscriptQuantification> listOfTranscriptQuantifications, String outmatrix) throws IOException {
		Path p = Paths.get(outmatrix);
		BufferedWriter bf = Files.newBufferedWriter(p,Charset.defaultCharset());
		
		int index = 0;
		for(TranscriptQuantification tq : listOfTranscriptQuantifications) {
			if(index == 0) {
				bf.append(tq.printHeader());
				bf.append("\n");
			}
			index++;
			tq.appendTo(bf);
			bf.append("\n");
		}
		
		bf.close();
	}
}
