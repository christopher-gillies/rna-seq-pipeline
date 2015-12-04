package org.kidneyomics.gtf;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class TranscriptLengthQuantifier {

	Map<String,List<Feature>> transcriptMap;
	String type = "exon";
	
	private TranscriptLengthQuantifier() {
		this.transcriptMap = new HashMap<String,List<Feature>>();
	}
	
	private TranscriptLengthQuantifier(String type) {
		this.transcriptMap = new HashMap<String,List<Feature>>();
		this.type = type;
	}
	
	public static TranscriptLengthQuantifier getTranscriptLengthQuantifier(String type) {
		return new TranscriptLengthQuantifier(type);
	}
	
	public static TranscriptLengthQuantifier getTranscriptLengthQuantifier() {
		return new TranscriptLengthQuantifier();
	}
	
	
	public void addExon(Feature feature) {
		String type = feature.type();
		if(type.equals(this.type)) {
			
			String transcriptId = feature.getAttribute("transcript_id");
			if(transcriptId == null) {
				throw new IllegalArgumentException("Feature has not transcript id");
			}
			
			if(!transcriptMap.containsKey(transcriptId)) {
				transcriptMap.put(transcriptId, new LinkedList<Feature>());
			} 
			
			transcriptMap.get(transcriptId).add(feature);
			
			
			
		} else {
			throw new IllegalArgumentException("TranscriptLengthQuantifier set up for " + type);
		}
	}
	
	/**
	 * 
	 * @param transcriptId
	 * @return the length of the transcript in terms the base pairs in all its exons. returns -1 if the transcript is not found. returns -2 if transcript has an invalid structure
	 */
	public int length(String transcriptId) {
		
		if(!transcriptMap.containsKey(transcriptId)) {
			return -1;
		}
		
		if(!validate(transcriptId)) {
			return -2;
		}
		
		List<Feature> features = transcriptMap.get(transcriptId);
		int length = 0;
		for(Feature f : features) {
			length += f.location().length();
		}
		return length;
	}
	
	/**
	 * 
	 * @param transcriptId
	 * @return true if exons do not overlap each other, false otherwise. false if no exons
	 */
	public boolean validate(String transcriptId) {
		List<Feature> features = transcriptMap.get(transcriptId);
		if(features.size() == 1) {
			return true;
		} else if(features.size() == 0) {
			return false;
		}
		
		Collections.sort(features, new FeatureComparator());
		Iterator<Feature> iter = features.iterator();
		
		Feature a = iter.next();
		while(iter.hasNext()) {
			Feature b = a;
			a = iter.next();
			if(a.location().intersection(b.location()) != null) {
				return false;
			}
		}
		return true;
	}
	
	
}
