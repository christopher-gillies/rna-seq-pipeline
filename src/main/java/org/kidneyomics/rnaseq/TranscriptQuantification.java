package org.kidneyomics.rnaseq;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.kidneyomics.gtf.FeatureComparator;
import org.kidneyomics.util.Chr2Int;
import org.springframework.util.StringUtils;

public class TranscriptQuantification implements Comparable<TranscriptQuantification>, Quantification {

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
	 * transcription_start_site
	 * start
	 * end
	 * length
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
	private int length = -1;
	private QuantificationType quantificationType = QuantificationType.TRANSCRIPT;
	private static FeatureComparator comparator = new FeatureComparator();
	
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
	
	public int getTranscriptionStartSite() {
		if(this.getStrand() == '+') {
			return this.feature.location().bioStart();
		} else {
			return this.feature.location().bioEnd();
		}
		
	}
	
	public int getStart() {
		return this.feature.location().bioStart();
	}
	
	public int getEnd() {
		return this.feature.location().bioEnd();
	}
	
	public int getLength() {
		return this.length;
	}
	
	public void setLength(int length) {
		this.length = length;
	}
	
	public char getStrand() {
		return this.feature.location().bioStrand();
	}



	public double getSampleExpression(String sample) {
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
		try {
			appendTo(sb);
		} catch(Exception e) {
			
		}
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
		
		sb.append("transcription_start_site");
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
	
	public void appendTo(Appendable appendable) throws IOException {
		
		appendable.append(getTranscriptId());
		appendable.append("\t");
		appendable.append(getGeneId());
		appendable.append("\t");
		appendable.append(getGeneName());
		appendable.append("\t");
		appendable.append(getGeneType());
		appendable.append("\t");
		appendable.append(getTranscriptType());
		appendable.append("\t");
		appendable.append(getChr());
		appendable.append("\t");
		appendable.append(Integer.toString(getTranscriptionStartSite()));
		appendable.append("\t");
		appendable.append(Integer.toString(getStart()));
		appendable.append("\t");
		appendable.append(Integer.toString(getEnd()));
		appendable.append("\t");
		appendable.append(Integer.toString(getLength()));
		appendable.append("\t");
		appendable.append(getStrand());
		appendable.append("\t");
		
		
		int index = 0;
		int last = sampleIds.size() - 1;
		
		for(String id : sampleIds) {
			appendable.append(Double.toString(getSampleExpression(id)));
			if(index != last) {
				appendable.append("\t");
			}
			index++;
		}
		
	}


	@Override
	public int compareTo(TranscriptQuantification o) {
		return comparator.compare(this.feature, o.feature);
	}
	
	public List<String> getSampleIds() {
		return this.sampleIds;
	}


	@Override
	public QuantificationType getQuantificationType() {
		return this.quantificationType;
	}
}
