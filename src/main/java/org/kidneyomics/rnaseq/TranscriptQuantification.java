package org.kidneyomics.rnaseq;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.springframework.util.StringUtils;

public class TranscriptQuantification {

	/*
	 * _______
	 * 
	 * COLUMNS
	 * _______
	 * transcript_id
	 * gene_id
	 * gene_symbol
	 * gene_type
	 * transcript_type
	 * chr
	 * start
	 * end
	 * strand
	 * sample1
	 * sample2
	 * ...
	 * samplen
	 * 
	 * 
	 * 
	 * 
	 */
	private Feature feature;
	private Map<String,Double> expressionMap;
	private List<String> sampleIds;
	
	public TranscriptQuantification(Feature feature, List<String> sampleIds) {
		this.feature = feature;
		this.sampleIds = sampleIds;
		this.expressionMap = new HashMap<String,Double>(sampleIds.size());
	}


	public Feature getFeature() {
		return feature;
	}
	
	public String getTranscriptId() {
		return this.feature.getAttribute("transcript_id");
	}
	
	public String getGeneId() {
		return this.feature.getAttribute("gene_id");
	}
	
	public String getGeneName() {
		return this.feature.getAttribute("gene_name");
	}
	
	public String getGeneType() {
		return this.feature.getAttribute("gene_type");
	}
	
	public String getTranscriptType() {
		return this.feature.getAttribute("transcript_type");
	}
	
	public String getChr() {
		return this.feature.seqname();
	}
	
	public int getStart() {
		return this.feature.location().bioStart();
	}
	
	public int getEnd() {
		return this.feature.location().bioEnd();
	}
	
	public int getLength() {
		return this.feature.location().length();
	}
	
	public char getStrand() {
		return this.feature.location().bioStrand();
	}



	public Double getSampleExpression(String sample) {
		return this.expressionMap.get(sample);
	}
	
	public void putSampleExpression(String sample, double expression) {
		this.expressionMap.put(sample, expression);
	}
	
	public void putSampleCount(String sample, Feature transcriptExpressionFeature) {
		this.expressionMap.put(sample, Double.parseDouble(transcriptExpressionFeature.getAttribute("reads")));
	}
	
	public void putSampleRPKM(String sample, Feature transcriptExpressionFeature) {
		this.expressionMap.put(sample, Double.parseDouble(transcriptExpressionFeature.getAttribute("RPKM")));
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		toString(sb);
		return sb.toString();
	}
	
	public String printHeader() {
		StringBuilder sb = new StringBuilder();
		sb.append("transcript_id");
		sb.append("\t");
		
		sb.append("gene_id");
		sb.append("\t");
		
		sb.append("gene_name");
		sb.append("\t");
		
		sb.append("gene_type");
		sb.append("\t");
		
		sb.append("transcript_type");
		sb.append("\t");
		
		sb.append("chr");
		sb.append("\t");
		
		sb.append("start");
		sb.append("\t");
		
		sb.append("end");
		sb.append("\t");
		
		sb.append("length");
		sb.append("\t");
		
		sb.append("strand");
		sb.append("\t");
		
		sb.append( StringUtils.collectionToDelimitedString(sampleIds, "\t")  );
	
		return sb.toString();
	}
	
	public void toString(StringBuilder sb) {
		
		sb.append(getTranscriptId());
		sb.append("\t");
		sb.append(getGeneId());
		sb.append("\t");
		sb.append(getGeneName());
		sb.append("\t");
		sb.append(getGeneType());
		sb.append("\t");
		sb.append(getTranscriptType());
		sb.append("\t");
		sb.append(getChr());
		sb.append("\t");
		sb.append(getStart());
		sb.append("\t");
		sb.append(getEnd());
		sb.append("\t");
		sb.append(getLength());
		sb.append("\t");
		sb.append(getStrand());
		sb.append("\t");
		
		
		int index = 0;
		int last = sampleIds.size() - 1;
		
		for(String id : sampleIds) {
			sb.append(getSampleExpression(id));
			if(index != last) {
				sb.append("\t");
			}
			index++;
		}
		
	}
	
	
}
