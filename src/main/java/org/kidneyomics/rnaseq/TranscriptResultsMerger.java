package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GeneLengthQuantifier;
import org.kidneyomics.gtf.TranscriptLengthQuantifier;
import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.slf4j.Logger;

abstract class TranscriptResultsMerger implements ApplicationCommand {

	protected ApplicationOptions applicationOptions;
	
	protected Logger logger;
	
	
	TranscriptResultsMerger(LoggerService loggerService, ApplicationOptions applicationOptions) {
		this.applicationOptions = applicationOptions;
		this.logger = loggerService.getLogger(this);
	}
	
	
	@Override
	public void doWork() throws Exception {
		Mode mode = applicationOptions.getMode();
		switch(mode) {
		case GENE_EXPRESSION_MATRIX: {
        	logger.info("Creating gene expression matrix");
        	this.writeGeneMatrix();
        	break;
        }
        case TRANSCRIPT_EXPRESSION_MATRIX: {
        	logger.info("Creating transcript expression matrix");
        	this.writeTranscriptMatrix();
        	break;
        }
        case TRANSCRIPT_RATIO_MATRIX:
        	logger.info("Creating transcript ratio matrix");
        	this.writeTranscriptRatioMatrix();
        	break;
        	default:
        		throw new IllegalArgumentException(mode + " not supported");
		}
		
	}

	public void writeTranscriptMatrix() throws Exception {
		String gtfList = applicationOptions.getFileIn();
		String annotation = applicationOptions.getGtf();
		String outmatrix = applicationOptions.getFileOut();
		boolean outCounts = applicationOptions.isOutCounts();
		Collection<TranscriptQuantification> listOfTranscriptQuantifications = getTranscriptQuantifications(gtfList, annotation, outCounts).getTranscriptQuantifications();
		logger.info("Writing expression matrix");
		QuantificationUtils.writeQuantificationMatrix(listOfTranscriptQuantifications,outmatrix);
		
	}
	
	public void writeGeneMatrix() throws Exception {
		String gtfList = applicationOptions.getFileIn();
		String annotation = applicationOptions.getGtf();
		String outmatrix = applicationOptions.getFileOut();
		boolean outCounts = applicationOptions.isOutCounts();
		TranscriptQuantificationResult transcriptQuantificationResult = getTranscriptQuantifications(gtfList, annotation, outCounts);
		List<GeneQuantification> listOfGeneQuantifications = getGeneQuantifications(transcriptQuantificationResult);
		QuantificationUtils.writeQuantificationMatrix(listOfGeneQuantifications,outmatrix);	
	}
	
	public void writeTranscriptRatioMatrix() throws Exception {
		String outmatrix = applicationOptions.getFileOut();
		List<TranscriptQuantification> tqs = getTranscriptRatios().transcriptRatios;
		QuantificationUtils.writeQuantificationMatrix(tqs,outmatrix);	
	}

	
	List<GeneQuantification> getGeneQuantifications(TranscriptQuantificationResult transcriptQuantificationResult) {
		
		List<TranscriptQuantification> listOfTranscriptQuantifications = transcriptQuantificationResult.getTranscriptQuantifications();
		List<GeneQuantification> list = new ArrayList<GeneQuantification>(transcriptQuantificationResult.getTranscriptQuantifications().size());
		HashMap<String,List<TranscriptQuantification>> quantificationsByGene = new HashMap<String,List<TranscriptQuantification>>();
		
		//calculates the length of a gene
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
			
			//transcript quantifications for this gene
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
		logger.info("Sorting results..");
		Collections.sort(list);
		return list;
	}
	
	TranscriptRatioResult getTranscriptRatios() throws Exception {
		logger.info("Getting transcript ratios");
		String gtfList = applicationOptions.getFileIn();
		String annotation = applicationOptions.getGtf();
		boolean outCounts = applicationOptions.isOutCounts();
		TranscriptQuantificationResult transcriptQuantificationResult = getTranscriptQuantifications(gtfList, annotation, outCounts);
		List<GeneQuantification> listOfGeneQuantifications = getGeneQuantifications(transcriptQuantificationResult);
		
		TranscriptRatioQuantifier trq = new TranscriptRatioQuantifier(listOfGeneQuantifications, transcriptQuantificationResult.sampleIds);
		List<TranscriptQuantification> tqs = trq.updateTranscriptRatios();
		
		TranscriptRatioResult trr = new TranscriptRatioResult(transcriptQuantificationResult, tqs);
		return trr;
	}
	
	/**
	 * @param fileIn should be a list of samples and their corresponding expression tables
	 * @return list of sample ids for some quantifications
	 */
	abstract List<String> getReadFileAndGetSampleIds(String fileIn) throws IOException;
	
	/**
	 * add transcript quantifications for each sample from file read in getReadFileAndGetSampleIds
	 * @param transcripts add transcript quantifications for each sample
	 * @param outCounts set if this is a count of an rpkm
	 * @throws Exception
	 */
	abstract void addTranscriptQuantifications(Map<String,TranscriptQuantification> transcripts, boolean outCounts) throws Exception;
	
	TranscriptQuantificationResult getTranscriptQuantifications(String fileIn, String annotation, boolean outCounts) throws Exception {
		
		/*
		 * Read sample information
		 */
		List<String> sampleIds = getReadFileAndGetSampleIds(fileIn);
		Collections.sort(sampleIds);
		
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
		
		
		logger.info("Obtaining transcript quantification results");
		addTranscriptQuantifications(transcripts, outCounts);
		
		
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
		return new TranscriptQuantificationResult(list,transcriptLengthQuantifier,geneLengthQuantifier,sampleIds);
	}
	
	
	protected class TranscriptQuantificationResult {
		List<TranscriptQuantification> transcriptQuantifications;
		TranscriptLengthQuantifier transcriptLengthQuantifier;
		GeneLengthQuantifier geneLengthQuantifier;
		List<String> sampleIds;
		
		TranscriptQuantificationResult(List<TranscriptQuantification> transcriptQuantifications,
				TranscriptLengthQuantifier transcriptLengthQuantifier,
				GeneLengthQuantifier geneLengthQuantifier, List<String> sampleIds) {
			this.transcriptLengthQuantifier = transcriptLengthQuantifier;
			this.transcriptQuantifications = transcriptQuantifications;
			this.geneLengthQuantifier = geneLengthQuantifier;
			this.sampleIds = sampleIds;
		}

		List<TranscriptQuantification> getTranscriptQuantifications() {
			return transcriptQuantifications;
		}

		void setTranscriptQuantifications(List<TranscriptQuantification> transcriptQuantifications) {
			this.transcriptQuantifications = transcriptQuantifications;
		}

		TranscriptLengthQuantifier getTranscriptLengthQuantifier() {
			return transcriptLengthQuantifier;
		}

		void setTranscriptLengthQuantifier(TranscriptLengthQuantifier transcriptLengthQuantifier) {
			this.transcriptLengthQuantifier = transcriptLengthQuantifier;
		}

		GeneLengthQuantifier getGeneLengthQuantifier() {
			return geneLengthQuantifier;
		}

		void setGeneLengthQuantifier(GeneLengthQuantifier geneLengthQuantifier) {
			this.geneLengthQuantifier = geneLengthQuantifier;
		}
		
		
	}
	
	protected class TranscriptRatioResult {
		TranscriptQuantificationResult transcriptQuantificationResult;
		List<TranscriptQuantification> transcriptRatios;
		
		TranscriptRatioResult(TranscriptQuantificationResult transcriptQuantificationResult,List<TranscriptQuantification> transcriptRatios) {
			this.transcriptQuantificationResult = transcriptQuantificationResult;
			this.transcriptRatios = transcriptRatios;
		}
	}

}
