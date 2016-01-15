package org.kidneyomics.rnaseq.stats;

import java.util.HashMap;
import java.util.Map;

class ModeCalculator {
	private HashMap<Integer,Integer> map;
	
	ModeCalculator() {
		this.map = new HashMap<>();
	}
	
	void add(int i) {
		if(map.containsKey(i)) {
			int newVal = map.get(i) + 1;
			map.put(i, newVal);
		} else {
			map.put(i, 1);
		}
	}
	
	int getMode() {
		int max = Integer.MIN_VALUE;
		int maxLoc = Integer.MIN_VALUE;
		for(Map.Entry<Integer, Integer> entry : map.entrySet()) {
			if(entry.getValue() > max) {
				max = entry.getValue();
				maxLoc = entry.getKey();
			}
		}
		return maxLoc;
	}
}
