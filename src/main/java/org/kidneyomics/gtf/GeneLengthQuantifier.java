package org.kidneyomics.gtf;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class GeneLengthQuantifier {

	Map<String,List<Feature>> geneMap;
	String type = "exon";
	
	private GeneLengthQuantifier() {
		this.geneMap = new HashMap<String,List<Feature>>();
	}
	
	private GeneLengthQuantifier(String type) {
		this.geneMap = new HashMap<String,List<Feature>>();
		this.type = type;
	}
	
	public static GeneLengthQuantifier getGeneLengthQuantifier(String type) {
		return new GeneLengthQuantifier(type);
	}
	
	public static GeneLengthQuantifier getGeneLengthQuantifier() {
		return new GeneLengthQuantifier();
	}
	
	
	public void addExon(Feature feature) {
		String type = feature.type();
		if(type.equals(this.type)) {
			
			String geneId = feature.getAttribute("gene_id");
			if(geneId == null) {
				throw new IllegalArgumentException("Feature has not gene id");
			}
			
			if(!geneMap.containsKey(geneId)) {
				geneMap.put(geneId, new LinkedList<Feature>());
			} 
			
			geneMap.get(geneId).add(feature);
			
			
			
		} else {
			throw new IllegalArgumentException("GeneLengthQuantifier set up for " + type);
		}
	}
	
	/**
	 * 
	 * @param geneId
	 * @return the length of the gene in terms the base pairs in all its exons (overlapping exons unioned together). returns -1 if the gene is not found. returns -2 if gene has an invalid structure
	 */
	public int length(String geneId) {
		
		if(!geneMap.containsKey(geneId)) {
			return -1;
		}
		
		if(!validate(geneId)) {
			return -2;
		}
		
		List<Feature> features = geneMap.get(geneId);
		List<Feature> merged = FeatureMerger.mergeOverlappingFeatures(features);
		
		int length = 0;
		for(Feature f : merged) {
			length += f.location().length();
		}
		return length;
	}
	
	/**
	 * 
	 * @param geneId
	 * @return true if exons do not overlap each other, false otherwise. false if no exons
	 */
	public boolean validate(String geneId) {
		List<Feature> features = geneMap.get(geneId);
		if(features.size() == 1) {
			return true;
		} else if(features.size() == 0) {
			return false;
		}
		
		Feature first = features.get(0);
		for(Feature feature : features) {
			if(first.location().bioStrand() != feature.location().bioStrand()) {
				return false;
			}
		}
		
		return true;
	}
	
	
}
