package org.kidneyomics.rnaseq.stats;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMRecord;

/**
 * TODO: This class is not complete; it may have to be rewritten to only look at k-mers found within the sequences
 * @author cgillies
 *
 */
class KmerStatistic extends AbstractReadPairStatistic {

	private final Map<String,Integer> kmerCount;
	private final int k;
	
	public KmerStatistic(int k) {
		super();		
		this.k = k;
		this.kmerCount = new HashMap<>();		
	}
	
	@Override
	protected void addRecord(SAMRecord record) {
		
		//get read string
		String bases = record.getReadString();
		
		//check all kmers of length k
		KmerIterator iter = new KmerIterator(k, bases);
		while(iter.hasNext()) {
			String kmer = iter.next();
			if(kmerCount.containsKey(kmer)) {
				int val = kmerCount.get(kmer);
				kmerCount.put(kmer, val + 1);
			} else {
				kmerCount.put(kmer, 1);
			}
		}

	}

	@Override
	public List<Double> getStatistic() {
		String[] fields = getFields();
		ArrayList<Double> results = new ArrayList<>(fields.length);
		for(int i = 0; i < fields.length; i++) {
			double val = kmerCount.get(fields[i]);
			results.add(val);
		}
		return results;
	}

	@Override
	protected String[] getFields() {
		String[] fields = new String[kmerCount.keySet().size()];
		kmerCount.keySet().toArray(fields);
		return fields;
	}
	
}
