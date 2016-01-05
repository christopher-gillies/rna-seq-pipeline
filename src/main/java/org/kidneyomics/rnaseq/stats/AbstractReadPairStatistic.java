package org.kidneyomics.rnaseq.stats;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.kidneyomics.rnaseq.SAMRecordPair;

import htsjdk.samtools.SAMRecord;

abstract class AbstractReadPairStatistic implements ReadPairStatistic {

	//the default approach is to apply the addRecord function to each record separately
	@Override
	public void addReadPair(SAMRecordPair pair) {
		SAMRecord pair1 = pair.getMate1();
		SAMRecord pair2 = pair.getMate2();
		addRecord(pair1);
		addRecord(pair2);
	}
	
	/**
	 * this should be implemented by subclass it specifies how to handle each record in the record pair
	 * @param record
	 */
	protected abstract void addRecord(SAMRecord record);
	
	/**
	 * 
	 * @return an array containing the fields to be tab separated. THIS MUST BE IN THE SAME ORDER AS THE getStatistic() methods
	 */
	protected abstract String[] getFields();
	
	@Override
	public String getHeader() {
		return StringUtils.join(getFields(), '\t');
	}

	@Override
	public abstract List<Double> getStatistic();

	@Override
	public void appendStatistic(Appendable appendable) throws IOException {
		appendable.append(StringUtils.join(getStatistic(), '\t'));
	}

	@Override
	public void appendHeader(Appendable appendable) throws IOException {
		appendable.append(getHeader());
	}
	
	/**
	 * 
	 */
	@Override
	public Map<String,Double> getStatisticAsMap() {
		
		String[] keys = getFields();
		List<Double> values = getStatistic();
		
		if(keys.length != values.size()) {
			throw new RuntimeException("The statistic has an unequal number of elements");
		}
		
		HashMap<String,Double> map = new HashMap<>();
		
		Iterator<Double> iter = values.iterator();
		int i = 0;
		
		while(iter.hasNext()) {
			Double value = iter.next();
			String key = keys[i];
			map.put(key, value);
			i++;
		}
		
		return map;
		
	}
}
