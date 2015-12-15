package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class GTFFeatureRenderer {

	public static String render(Feature feature) {
		StringBuilder sb = new StringBuilder();
		
		sb.append(feature.seqname());
		sb.append("\t");
		
		sb.append(feature.source());
		sb.append("\t");
		
		sb.append( feature.type() );
		sb.append("\t");
		
		sb.append( feature.location().bioStart()  );
		sb.append("\t");
		
		sb.append( feature.location().bioEnd() );
		sb.append("\t");
		
		if(feature.score() != Double.MIN_VALUE && feature.score() != Double.MAX_VALUE) {
			sb.append(feature.score() );
		} else {
			sb.append(".");
		}
		sb.append("\t");
		
		sb.append( feature.location().bioStrand() );
		sb.append("\t");
		
		if(feature.frame() != Integer.MIN_VALUE && feature.frame() != Integer.MAX_VALUE) {
			sb.append( feature.frame() );
		} else {
			sb.append(".");
		}
		
		
		sb.append("\t");
		
		sb.append( feature.attributes() );
		
		return sb.toString();
		
	}
	
	
}
