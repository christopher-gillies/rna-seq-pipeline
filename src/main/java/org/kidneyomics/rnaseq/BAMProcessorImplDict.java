package org.kidneyomics.rnaseq;

public class BAMProcessorImplDict implements BAMProcessor {

	@Override
	public SAMRecordPair getNextReadPair() {
		//insert first read into dictionary by queryname
		//insert second read into dictionary
		//check if the dictionary length for that entry has both pairs
		//if it is does return the read pair
		//otherwise continue reading
		//this way just return pairs as they are completed
		//should be MUCH faster
		return null;
	}

}
