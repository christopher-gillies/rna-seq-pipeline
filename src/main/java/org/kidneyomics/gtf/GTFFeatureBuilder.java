package org.kidneyomics.gtf;


import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.Location;


public class GTFFeatureBuilder {
	
	/**
	 * 
	 * @param line
	 * @return a biojava Feature object representing the line
	 */
	public static Feature createFromLine(String line) {
		return createFromLine(line,false);
	}
	
	/**
	 * 
	 * @param line
	 * @param removeVersion - true or false to remove the version number from a transcript/gene/exon identifier
	 * @return a biojava Feature object representing the line
	 */
	public static Feature createFromLine(String line, boolean removeVersion) {
		String[] cols = line.split("\t");
		if(cols.length != 9) {
			throw new IllegalArgumentException("Incorrect number of columns. Found " + cols.length);
		}
		String seqname = cols[0];
		String source = cols[1];
		String type = cols[2];
		char strand = cols[6].charAt(0);
		Location location = Location.fromBio(Integer.parseInt(cols[3]), Integer.parseInt(cols[4]), strand);
		
		Double score = null;
		if(cols[5].equalsIgnoreCase(".")) {
			score = Double.MIN_VALUE;
		} else {
			score = Double.parseDouble(cols[5]);
		}
		
		Integer frame = null;
		
		if(cols[7].equals(".")) {
			frame = Integer.MIN_VALUE;
		} else {
			frame = Integer.parseInt(cols[7]);
		}
		
		String attributes = cols[8];
		//return new GTFFeature(seqname, source, type, location, score, strand, frame, attributes);
		
		if(removeVersion) {
			attributes = VersionTrimmer.trim(attributes);
		}
		
		return new Feature(seqname, source, type, location, score, frame, attributes);
	}
	
	public static FeatureList createFeatureList(List<String> lines) {
		FeatureList list = new FeatureList();
		for(String line : lines) {
			list.add(createFromLine(line));
		}
		return list;
	}
	
	public static Feature addAttributesToFeature(Feature feature, Map<String,String> items) {
		StringBuilder sb = new StringBuilder();
		sb.append(feature.attributes());
		for(Map.Entry<String, String> item : items.entrySet()) {
			sb.append(" ");
			sb.append(item.getKey());
			sb.append(" \"");
			sb.append(item.getValue());
			sb.append("\"");
			sb.append(";");
		}
		String newAtt = sb.toString();
		
		return new Feature( feature.seqname(), feature.source(), feature.type(), feature.location(), feature.score(), feature.frame(), newAtt  );
	}
}
