package org.kidneyomics.rnaseq;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.FeatureComparator;
import org.kidneyomics.gtf.GTF_FORMAT;
import org.springframework.util.StringUtils;

class ExonOrGeneQuantificationResult implements MutableQuantification {

	/*
	 * _______
	 * 
	 * COLUMNS
	 * _______
	 * [Optional exon_id]
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
	
	private ExonOrGeneQuantificationResult(Feature feature, List<String> sampleIds, QuantificationType type) {
		this.feature = feature;
		this.sampleIds = sampleIds;
		this.expressionMap = new HashMap<String,Double>(sampleIds.size());
		this.quantificationType = type;
		this.format = GTF_FORMAT.getFormat(feature);
		
		String fType = feature.type();
		//Set the length
		switch(this.quantificationType) {
		case GENE:
		case GENE_COUNTS:
		case GENE_RPKM:	
			if(!fType.equals("gene")) {
				throw new IllegalArgumentException("Quantification type does not match feature type");
			}
			if(this.feature.location().bioStrand() == '+') {
				this.tss = feature.location().bioStart();
			} else {
				this.tss = feature.location().bioEnd();
			}
			this.length = Integer.parseInt(feature.getAttribute("length"));
			break;
		case EXON:
		case EXON_COUNTS:
		case EXON_RPKM:
			if(!fType.equals("exon")) {
				throw new IllegalArgumentException("Quantification type does not match feature type");
			}
			this.length = feature.location().length();
			this.tss = Integer.parseInt(feature.getAttribute("tss"));
			break;
			default:
				this.length = -1;
				this.tss = -1;
		}
	}

	public static ExonOrGeneQuantificationResult getExonCountsQuantification(Feature feature, List<String> sampleIds) {
		return new ExonOrGeneQuantificationResult(feature,sampleIds,QuantificationType.EXON_COUNTS);
	}
	
	public static ExonOrGeneQuantificationResult getExonRPKMQuantification(Feature feature, List<String> sampleIds) {
		return new ExonOrGeneQuantificationResult(feature,sampleIds,QuantificationType.EXON_RPKM);
	}
	
	public static ExonOrGeneQuantificationResult getGeneCountsQuantification(Feature feature, List<String> sampleIds) {
		return new ExonOrGeneQuantificationResult(feature,sampleIds,QuantificationType.GENE_COUNTS);
	}
	
	public static ExonOrGeneQuantificationResult getGeneRPKMQuantification(Feature feature, List<String> sampleIds) {
		return new ExonOrGeneQuantificationResult(feature,sampleIds,QuantificationType.GENE_RPKM);
	}
	
	private final Feature feature;
	private final Map<String,Double> expressionMap;
	private final List<String> sampleIds;
	private int length;
	private final int tss;
	private final GTF_FORMAT format;
	private final QuantificationType quantificationType;
	private static final FeatureComparator comparator = new FeatureComparator();
	
	
	@Override
	public String getGeneId() {
		return this.feature.getAttribute("gene_id");
	}
	
	@Override
	public String getGeneName() {
		return this.feature.getAttribute("gene_name");
	}
	
	@Override
	public String getGeneType() {
		if(this.format == GTF_FORMAT.GENCODE) {
			return this.feature.getAttribute("gene_type");
		} else {
			return this.feature.getAttribute("gene_biotype");
		}
	}
	
	@Override
	public String getChr() {
		return this.feature.seqname();
	}
	
	@Override
	public int getTranscriptionStartSite() {
		return this.tss;
	}
	
	@Override
	public int getStart() {
		return this.feature.location().bioStart();
	}
	
	@Override
	public int getEnd() {
		return this.feature.location().bioEnd();
	}
	
	@Override
	public int getLength() {
		return this.length;
	}
	
	@Override
	public void setLength(int length) {
		this.length = length;
	}
	
	@Override
	public char getStrand() {
		return this.feature.location().bioStrand();
	}

	@Override
	public double getSampleExpression(String sample) {
		
		Double result = this.expressionMap.get(sample);
		if(result == null) {
			return 0.0;
		} else {
			return result;
		}
	}
	
	@Override
	public void putSampleExpression(String sample, double expression) {
		this.expressionMap.put(sample, expression);
	}
	
	@Override
	public void putSampleCount(String sample, Feature transcriptExpressionFeature) {
		if(this.getQuantificationType() != QuantificationType.EXON_COUNTS && this.getQuantificationType() != QuantificationType.GENE_COUNTS) {
			throw new IllegalAccessError("Only can add count values for this feature");
		}
		this.expressionMap.put(sample, Double.parseDouble(transcriptExpressionFeature.getAttribute("reads")));
	}
	
	@Override
	public void putSampleRPKM(String sample, Feature transcriptExpressionFeature) {
		if(this.getQuantificationType() != QuantificationType.EXON_RPKM && this.getQuantificationType() != QuantificationType.GENE_RPKM) {
			throw new IllegalAccessError("Only can add RPKM values for this feature");
		}
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
	
	@Override
	public String printHeader() {
		StringBuilder sb = new StringBuilder();
		if(this.quantificationType == QuantificationType.EXON_COUNTS || this.quantificationType == QuantificationType.EXON_RPKM) {
			sb.append("exon_id");
			sb.append("\t");
		}
		
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
	
	@Override
	public void appendTo(Appendable appendable) throws IOException {
		
		if(this.quantificationType == QuantificationType.EXON_COUNTS || this.quantificationType == QuantificationType.EXON_RPKM) {
			appendable.append(this.feature.getAttribute("id"));
			appendable.append("\t");
		}
		
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


	GTF_FORMAT getFORMAT() {
		return this.format;
	}
	
	
	public List<String> getSampleIds() {
		return this.sampleIds;
	}


	@Override
	public QuantificationType getQuantificationType() {
		return this.quantificationType;
	}

	@Override
	public int compareTo(Quantification o) {
		return comparator.compare(this.feature, o.getFeature());
	}

	@Override
	public String getId() {
		if(this.quantificationType == QuantificationType.EXON_COUNTS || this.quantificationType == QuantificationType.EXON_RPKM) {
			return this.feature.getAttribute("id");
		} else {
			return this.getGeneId();
		}
	}

	@Override
	public Feature getFeature() {
		return this.feature;
	}

}
