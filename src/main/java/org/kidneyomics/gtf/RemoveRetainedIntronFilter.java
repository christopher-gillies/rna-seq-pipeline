package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class RemoveRetainedIntronFilter implements GTFFeatureFilter {

	private final String retained_intron_string = "retained_intron";
	@Override
	public boolean filter(Feature feature) {

		String retained = "";
		if(feature.hasAttribute("transcript_type")) {
			retained = feature.getAttribute("transcript_type");
		}


		
		return !(retained.equals(retained_intron_string) || feature.source().equals(retained_intron_string));
	}

}
