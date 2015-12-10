package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class FeatureCount {
	
	private Feature feature;
	private double count = 0.0;
	private String id = null;
	
	public FeatureCount(Feature feature) {
		if(feature == null) {
			throw new IllegalArgumentException("feature cannot be null");
		}
		this.feature = feature;
	}
	
	public Feature getFeature() {
		return feature;
	}
	public void setFeature(Feature feature) {
		this.feature = feature;
	}
	
	public double getCount() {
		return count;
	}
	
	public void addToCount(double val) {
		assert(val >= 0);
		this.count += val;
	}
	
	public double getRPKM(double numberOfReads) {
		double length = feature.location().length();
		return ((double) count) / length / numberOfReads * Math.pow(10, 9);
	}
	
	public String getId() {
		if(id == null) {
			this.id = FeatureCount.featureIdMaker(this.feature);
			return this.id;
		} else {
			return this.id;
		}
	}
	
	public static String featureIdMaker(Feature feature) {
		String geneId = feature.getAttribute("gene_id");
		int start = feature.location().bioStart();
		int end = feature.location().bioEnd();
		
		StringBuilder sb = new StringBuilder();
		sb.append(geneId);
		sb.append("_");
		sb.append(feature.seqname());
		sb.append("_");
		sb.append(start);
		sb.append("_");
		sb.append(end);
		return sb.toString();
	}
	
	
}
