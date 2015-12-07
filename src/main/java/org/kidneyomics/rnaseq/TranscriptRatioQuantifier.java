package org.kidneyomics.rnaseq;

import java.util.LinkedList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

class TranscriptRatioQuantifier {

	private List<GeneQuantification> geneQuantifications;
	private List<String> sampleIds;
	
	private Logger logger = LoggerFactory.getLogger(TranscriptRatioQuantifier.class);
	
	TranscriptRatioQuantifier(List<GeneQuantification> geneQuantifications, List<String> sampleIds) {
		this.geneQuantifications = geneQuantifications;
		this.sampleIds = sampleIds;
	}
	
	List<TranscriptQuantification> updateTranscriptRatios() {
		
		List<TranscriptQuantification> transcriptQuantifications = new LinkedList<TranscriptQuantification>();
		
		for(GeneQuantification gq : geneQuantifications) {
			
			//logger.info("GENE: " +gq.getGeneId());
			/*
			 * Update expression for each sample for each transcript in gene
			 */
			List<TranscriptQuantification> tqs = gq.getTranscriptQuantifications();
			for(String id : sampleIds) {
				double totalExpression = gq.getSampleExpression(id);
				//logger.info("SAMPLE: " + id);
				//logger.info("GENE EXPRESSION: " + totalExpression);
				if(totalExpression != 0.0) {
					for(TranscriptQuantification tq : tqs) {
						double transcriptExpression = tq.getSampleExpression(id);
						//logger.info("SAMPLE: " + id);
						//logger.info("TRANSCRIPT EXPRESSION: " + transcriptExpression);
						//logger.info("TRANSCRIPT RATIO: " + transcriptExpression / totalExpression);
						tq.putSampleExpression(id, transcriptExpression / totalExpression);
					}
				}
			}
			
			//Add transcripts to transcriptQuantifications
			transcriptQuantifications.addAll(tqs);
		}
		
		return transcriptQuantifications;
	}
	
}
