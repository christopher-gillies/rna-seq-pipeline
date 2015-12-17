package org.kidneyomics.rnaseq;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.kidneyomics.gtf.FeatureComparator;
import org.springframework.util.StringUtils;

public class GeneQuantification implements Quantification {

	
	/*
	 * _______
	 * 
	 * COLUMNS
	 * _______
	 * gene_id
	 * gene_symbol
	 * gene_type
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
	
	private String geneId;
	private String geneName;
	private String geneType;
	private String chr;
	private int start;
	private int end;
	private int length;
	private char strand;
	private List<TranscriptQuantification> transcriptQuantifications;
	private Map<String,Double> expressionMap;
	private List<String> sampleIds;
	private static FeatureComparator comparator = new FeatureComparator();
	private Feature feature;
	private QuantificationType quantificationType = QuantificationType.GENE;
	
	public GeneQuantification(List<TranscriptQuantification> transcriptQuantifications) {
		if(transcriptQuantifications == null || transcriptQuantifications.size() == 0) {
			throw new IllegalArgumentException("A gene must have at least one transcript");
		}
		
		this.length = -1;
		
		this.geneId = transcriptQuantifications.get(0).getGeneId();
		this.geneName = transcriptQuantifications.get(0).getGeneName();
		this.geneType = transcriptQuantifications.get(0).getGeneType();
		this.chr = transcriptQuantifications.get(0).getChr();
		this.start = transcriptQuantifications.get(0).getStart();
		this.end = transcriptQuantifications.get(0).getEnd();
		this.strand = transcriptQuantifications.get(0).getStrand();
		this.sampleIds = transcriptQuantifications.get(0).getSampleIds();
		
		
		/*
		 * Make sure all transcripts are from the same gene
		 */
		for(TranscriptQuantification tq : transcriptQuantifications) {
			if(!tq.getGeneId().equals(this.geneId)) {
				throw new IllegalArgumentException("Transcripts must all be from the same gene");
			}
		}
		
		/*
		 * Find lowest start position and find highest end position
		 */
		for(TranscriptQuantification tq : transcriptQuantifications) {
			if(tq.getStart() < this.start) {
				this.start = tq.getStart();
			}
				
			if(tq.getEnd() > this.end) {
				this.end = tq.getEnd();
			}
		}
		
		/*
		 * Make sure all transcripts have the same samples
		 */
		for(TranscriptQuantification tq : transcriptQuantifications) {
			List<String> sampleIdsForTq = tq.getSampleIds();
			if(!sampleIdsForTq.containsAll(this.sampleIds)) {
				throw new IllegalArgumentException("Transcripts must all have the same sample ids");
			}
		}
		
		this.transcriptQuantifications = transcriptQuantifications;
		
		/*
		 * set feature
		 */
		Location location = Location.fromBio(this.start, this.end, this.strand);
		this.feature = new Feature(this.chr, "", "gene", location, 0.0, 0, "");
	}

	/**
	 * 
	 * @return the start position if + strand and the end position if - strand
	 */
	public int getTranscriptionStartSite() {
		if(this.getStrand() == '+') {
			return this.start;
		} else {
			return this.end;
		}
		
	}
	
	public String getGeneId() {
		return geneId;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getGeneType() {
		return geneType;
	}

	public String getChr() {
		return chr;
	}

	/**
	 * 
	 * @return the lowest start position among all transcripts
	 */
	public int getStart() {
		return start;
	}

	/**
	 * 
	 * @return the highest end position among all transcripts
	 */
	public int getEnd() {
		return end;
	}

	public char getStrand() {
		return strand;
	}


	/**
	 * 
	 * @return the length of the gene in terms of base pairs in exons. Overlapping exons are merged into super exons.
	 */
	public int getLength() {
		return length;
	}


	public void setLength(int length) {
		this.length = length;
	}
	
	protected void makeExpressionMap() {
		this.expressionMap = new HashMap<String,Double>();
		
		/*
		 * initialize all samples to 0
		 */
		for(String id : this.sampleIds) {
			this.expressionMap.put(id, 0.0);
		}
		
		/*
		 * compute the sum of expression across all transcripts for each sample individually 
		 */
		for(TranscriptQuantification tq : this.transcriptQuantifications) {
			for(String id : this.sampleIds) {
				double expression = tq.getSampleExpression(id);
				double newExpression = expression + this.expressionMap.get(id);
				//set expression
				this.expressionMap.put(id, newExpression);
			}
		}
	}
	
	/**
	 * 
	 * @return a map of the expression as the sum of each transcript's expression per sample
	 */
	protected Map<String,Double> getExpressionMap() {
		if(this.expressionMap != null) {
			return this.expressionMap;
		} else {
			makeExpressionMap();
			return this.expressionMap;
		}
	}
	
	/**
	 * 
	 * @param id
	 * @return the gene expression for sample id
	 */
	public double getSampleExpression(String id) {
		if(this.expressionMap != null) {
			return this.expressionMap.get(id);
		} else {
			makeExpressionMap();
			return this.expressionMap.get(id);
		}
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

		sb.append("gene_id");
		sb.append("\t");
		
		sb.append("gene_name");
		sb.append("\t");
		
		sb.append("gene_type");
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
		
		appendable.append(getGeneId());
		appendable.append("\t");
		appendable.append(getGeneName());
		appendable.append("\t");
		appendable.append(getGeneType());
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
	public int compareTo(Quantification o) {
		return comparator.compare(this.feature, o.getFeature());
	}

	@Override
	public QuantificationType getQuantificationType() {
		return this.quantificationType;
	}
	
	List<TranscriptQuantification> getTranscriptQuantifications() {
		return this.transcriptQuantifications;
	}

	@Override
	public String getId() {
		return getGeneId();
	}

	@Override
	public Feature getFeature() {
		return this.feature;
	}
	
	
}
