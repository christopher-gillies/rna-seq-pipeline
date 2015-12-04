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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GeneLengthQuantifier;
import org.kidneyomics.gtf.TranscriptLengthQuantifier;
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
		Collection<TranscriptQuantification> listOfTranscriptQuantifications = getTranscriptQuantifications(gtfList, annotation, outCounts).getTranscriptQuantifications();
		writeQuantificationMatrix(listOfTranscriptQuantifications,outmatrix);
		
	}
	
	public void writeGeneMatrix() throws Exception {
		String gtfList = applicationOptions.getFileIn();
		String annotation = applicationOptions.getGtf();
		String outmatrix = applicationOptions.getFileOut();
		boolean outCounts = applicationOptions.isOutCounts();
		TranscriptQuantificationResult transcriptQuantificationResult = getTranscriptQuantifications(gtfList, annotation, outCounts);
		List<GeneQuantification> listOfGeneQuantifications = getGeneQuantifications(transcriptQuantificationResult);
		writeQuantificationMatrix(listOfGeneQuantifications,outmatrix);	
	}
	
	protected List<GeneQuantification> getGeneQuantifications(TranscriptQuantificationResult transcriptQuantificationResult) {
		
		List<TranscriptQuantification> listOfTranscriptQuantifications = transcriptQuantificationResult.getTranscriptQuantifications();
		List<GeneQuantification> list = new ArrayList<GeneQuantification>(transcriptQuantificationResult.getTranscriptQuantifications().size());
		HashMap<String,List<TranscriptQuantification>> quantificationsByGene = new HashMap<String,List<TranscriptQuantification>>();
		GeneLengthQuantifier glq = transcriptQuantificationResult.getGeneLengthQuantifier();
		
		/*
		 * Organize transcripts by geneid
		 */
		logger.info("Finding transcripts for each gene");
		for(TranscriptQuantification tq : listOfTranscriptQuantifications) {
			String geneId = tq.getGeneId();
			if(quantificationsByGene.containsKey(geneId)) {
				quantificationsByGene.get(geneId).add(tq);
			} else {
				LinkedList<TranscriptQuantification> listForGene = new LinkedList<TranscriptQuantification>();
				listForGene.add(tq);
				quantificationsByGene.put(geneId, listForGene);
			}
		}
		
		/*
		 * Create gene quantifications from list of transcripts for each gene
		 */
		logger.info("Calculating gene level quantification");
		for(Map.Entry<String, List<TranscriptQuantification>> entry : quantificationsByGene.entrySet()) {
			List<TranscriptQuantification> listForGene = entry.getValue();
			String geneId = entry.getKey();
			GeneQuantification gq = new GeneQuantification(listForGene);
			/*
			 * Calculate length
			 */
			int length = glq.length(geneId);
			gq.setLength(length);
			
			/*
			 * add to list
			 */
			list.add(gq);
		}
		
		/*
		 * Sort
		 */
		logger.info("Sorting results");
		Collections.sort(list);
		return list;
	}
	
	protected TranscriptQuantificationResult getTranscriptQuantifications(String gtfList, String annotation, boolean outCounts) throws Exception {
		
		/*
		 * Read sample information
		 */
		List<SampleGTF> sampleGtfs = sampleGtfReader.readFile(new File(gtfList));
		List<String> sampleIds = sampleGtfReader.getIds(sampleGtfs);
		
		TranscriptLengthQuantifier transcriptLengthQuantifier = TranscriptLengthQuantifier.getTranscriptLengthQuantifier("exon");
		GeneLengthQuantifier geneLengthQuantifier = GeneLengthQuantifier.getGeneLengthQuantifier("exon");
		
		/*
		 * Read gencode transcripts
		 */
		logger.info("Reading annotation from specified gtf file");
		HashMap<String,TranscriptQuantification> transcripts = new HashMap<String,TranscriptQuantification>();
		GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(new File(annotation));
		for(Feature feature : reader) {
			/*
			 * If it is a transcript..
			 */
			if(feature.type().equals("transcript")) {
				String transcriptId = feature.getAttribute("transcript_id");
				
				transcripts.put(transcriptId, new TranscriptQuantification(feature,sampleIds));
			} else if(feature.type().equals("exon")) {
				transcriptLengthQuantifier.addExon(feature);
				geneLengthQuantifier.addExon(feature);
			}
		}
		reader.close();
		logger.info("Finished reading annotation");
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
		logger.info("Finished reading transcript expression");
		

		/*
		 * Get collection
		 */
		Collection<TranscriptQuantification> listOfTranscriptQuantifications = transcripts.values();
		
		/*
		 * Update transcript lengths to be the sum of the base pairs in a transcript
		 */
		logger.info("Calculating transcript lengths");
		for(TranscriptQuantification tq : listOfTranscriptQuantifications) {
			String id = tq.getTranscriptId();
			int lengthOfTranscript = transcriptLengthQuantifier.length(id);
			tq.setLength(lengthOfTranscript);
		}
				
		/*
		 * sort
		 */
		logger.info("Sorting results...");
		List<TranscriptQuantification> list = new ArrayList<TranscriptQuantification>(listOfTranscriptQuantifications);
		Collections.sort(list);
		
		/*
		 * return results
		 */
		return new TranscriptQuantificationResult(list,transcriptLengthQuantifier,geneLengthQuantifier);
	}
	
	protected void writeQuantificationMatrix(Collection<? extends Quantification> listOfQuantifications, String outmatrix) throws IOException {
		logger.info("Writing expression matrix");
		Path p = Paths.get(outmatrix);
		BufferedWriter bf = Files.newBufferedWriter(p,Charset.defaultCharset());
		
		int index = 0;
		for(Quantification tq : listOfQuantifications) {
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
	
	protected class TranscriptQuantificationResult {
		List<TranscriptQuantification> transcriptQuantifications;
		TranscriptLengthQuantifier transcriptLengthQuantifier;
		GeneLengthQuantifier geneLengthQuantifier;
		
		protected TranscriptQuantificationResult(List<TranscriptQuantification> transcriptQuantifications,
				TranscriptLengthQuantifier transcriptLengthQuantifier,
				GeneLengthQuantifier geneLengthQuantifier) {
			this.transcriptLengthQuantifier = transcriptLengthQuantifier;
			this.transcriptQuantifications = transcriptQuantifications;
			this.geneLengthQuantifier = geneLengthQuantifier;
		}

		public List<TranscriptQuantification> getTranscriptQuantifications() {
			return transcriptQuantifications;
		}

		public void setTranscriptQuantifications(List<TranscriptQuantification> transcriptQuantifications) {
			this.transcriptQuantifications = transcriptQuantifications;
		}

		public TranscriptLengthQuantifier getTranscriptLengthQuantifier() {
			return transcriptLengthQuantifier;
		}

		public void setTranscriptLengthQuantifier(TranscriptLengthQuantifier transcriptLengthQuantifier) {
			this.transcriptLengthQuantifier = transcriptLengthQuantifier;
		}

		public GeneLengthQuantifier getGeneLengthQuantifier() {
			return geneLengthQuantifier;
		}

		public void setGeneLengthQuantifier(GeneLengthQuantifier geneLengthQuantifier) {
			this.geneLengthQuantifier = geneLengthQuantifier;
		}
		
		
	}
}
