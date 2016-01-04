package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

public class STARLogTest {

	@Test
	public void testRegexWithPercent() {
		String a = "Uniquely mapped reads % |	85.62%";
		String res = a.replaceAll(STARLog.regexWithPercent, "_");
		String exp = "Uniquely_mapped_reads___|_85.62_";
		assertEquals(exp,res);
	}
	
	@Test
	public void testRegexNoPercent() {
		String a = "Uniquely mapped reads % |	85.62%";
		String res = a.replaceAll(STARLog.regexNoPercent, "_");
		String exp = "Uniquely_mapped_reads_%_|_85.62%";
		assertEquals(exp,res);
	}
	
	@Test
	public void testCreateSTARLog() throws Exception {
		String in = "                           Started job on |	Nov 24 13:11:41\n" + 
				"                             Started mapping on |	Nov 24 13:12:15\n" + 
				"                                    Finished on |	Nov 24 13:18:34\n" + 
				"       Mapping speed, Million of reads per hour |	303.09\n" + 
				"\n" + 
				"                          Number of input reads |	31908652\n" + 
				"                      Average input read length |	78\n" + 
				"                                    UNIQUE READS:\n" + 
				"                   Uniquely mapped reads number |	27319614\n" + 
				"                        Uniquely mapped reads % |	85.62%\n" + 
				"                          Average mapped length |	77.67\n" + 
				"                       Number of splices: Total |	7769525\n" + 
				"            Number of splices: Annotated (sjdb) |	7757722\n" + 
				"                       Number of splices: GT/AG |	7616623\n" + 
				"                       Number of splices: GC/AG |	124419\n" + 
				"                       Number of splices: AT/AC |	19037\n" + 
				"               Number of splices: Non-canonical |	9446\n" + 
				"                      Mismatch rate per base, % |	0.22%\n" + 
				"                         Deletion rate per base |	0.00%\n" + 
				"                        Deletion average length |	1.43\n" + 
				"                        Insertion rate per base |	0.00%\n" + 
				"                       Insertion average length |	1.17\n" + 
				"                             MULTI-MAPPING READS:\n" + 
				"        Number of reads mapped to multiple loci |	2644721\n" + 
				"             % of reads mapped to multiple loci |	8.29%\n" + 
				"        Number of reads mapped to too many loci |	18984\n" + 
				"             % of reads mapped to too many loci |	0.06%\n" + 
				"                                  UNMAPPED READS:\n" + 
				"       % of reads unmapped: too many mismatches |	0.00%\n" + 
				"                 % of reads unmapped: too short |	5.96%\n" + 
				"                     % of reads unmapped: other |	0.07%";
		
		String [] linesArr = in.split("\n");
		List<String> lines = Arrays.asList(linesArr);
		STARLog log = STARLog.createLogFromLines("id",lines);
		String res = log.getValue("Uniquely_mapped_reads_number");
		assertEquals("27319614",res);
		
		String res2 = log.getValue("Uniquely_mapped_reads_%");
		assertEquals("85.62",res2);
		
		assertEquals(27,log.getNumberOfKeys());
		
	}

}
