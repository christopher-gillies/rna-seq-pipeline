package org.kidneyomics.rnaseq.stats;

import java.util.Iterator;

/*
 * Java Strings is really inefficient for this operation
 */
/**
 * 
 * @author cgillies
 * This class iterates through all kmers of length k for a string
 */
class KmerIterator implements Iterator<String> {

	private final int k;
	private final String bases;
	private final int length;
	private int currentIndex = 0;
	KmerIterator(final int k, final String bases) {
		this.k = k;
		this.bases = bases;
		this.length = bases.length();
		
		if(k <= 0) {
			throw new IllegalArgumentException("k must be greater than 0");
		}
		
		if(k > length) {
			throw new IllegalArgumentException("k must be at most the length of the bases string");
		}
	}
	
	
	
	@Override
	public boolean hasNext() {
		if(currentIndex < length - k + 1) {
			return true;
		} else {
			return false;
		}
	}

	//count save a little execution time by not rechecking hasNext() but the issue is we could overflow on currentIndex if it is call to many times. this seem quite unlikely though
	@Override
	public String next() {
		if(hasNext()) {
			//end position is exclusive
			String kmer = bases.substring(currentIndex, currentIndex + k);
			currentIndex++;
			return kmer;
		} else {
			return null;
		}
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}
	
	

}
