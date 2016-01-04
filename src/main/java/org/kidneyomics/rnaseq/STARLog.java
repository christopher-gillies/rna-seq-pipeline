package org.kidneyomics.rnaseq;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

class STARLog {

	private Map<String,String> map = new HashMap<>();
	private String id;
	
	/*
	 *                                 
	 *                           Started job on |	Nov 24 13:11:41
                             Started mapping on |	Nov 24 13:12:15
                                    Finished on |	Nov 24 13:18:34
       Mapping speed, Million of reads per hour |	303.09

                          Number of input reads |	31908652
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	27319614
                        Uniquely mapped reads % |	85.62%
                          Average mapped length |	77.67
                       Number of splices: Total |	7769525
            Number of splices: Annotated (sjdb) |	7757722
                       Number of splices: GT/AG |	7616623
                       Number of splices: GC/AG |	124419
                       Number of splices: AT/AC |	19037
               Number of splices: Non-canonical |	9446
                      Mismatch rate per base, % |	0.22%
                         Deletion rate per base |	0.00%
                        Deletion average length |	1.43
                        Insertion rate per base |	0.00%
                       Insertion average length |	1.17
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	2644721
             % of reads mapped to multiple loci |	8.29%
        Number of reads mapped to too many loci |	18984
             % of reads mapped to too many loci |	0.06%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	5.96%
                     % of reads unmapped: other |	0.07%
	 */
	
	static final String regexNoPercent = "[ ]|\t";
	static final String regexWithPercent = "[ %]|\t";
	
	private STARLog() {
		
	}
	
	static STARLog createLogFromLines(String id, List<String> lines) throws Exception {
		STARLog log = new STARLog();
		log.id = id;
		for(String line : lines) {
			if(line.contains("|")) {
				String[] cols = line.split("[|]");
				if(cols.length != 2) {
					//System.err.println(line);
					//for(int i = 0; i < cols.length; i++) {
					//	System.err.println(cols[i]);
					//}
					throw new Exception("STAR log has an unexpected format: " + "size: " + cols.length + " -- " + line);
				}
				log.map.put(cols[0].trim().replaceAll(regexNoPercent, "_"), cols[1].trim().replaceAll(regexWithPercent, ""));
			}
				
		}
		
		return log;
	}
	
	String getValue(String key) {
		String result = null;
		result = map.get(key);
		if(result == null) {
			result = "NA";
		}
		return result;
	}
	
	int getNumberOfKeys() {
		return this.map.size();
	}
	
	Set<String> getKeys() {
		return this.map.keySet();
	}
	
	String getId() {
		return this.id;
	}
}
