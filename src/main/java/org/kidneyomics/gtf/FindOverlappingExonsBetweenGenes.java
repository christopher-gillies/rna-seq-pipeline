package org.kidneyomics.gtf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
public class FindOverlappingExonsBetweenGenes {

	Logger logger = LoggerFactory.getLogger(FindOverlappingExonsBetweenGenes.class);
	/**
	 * 
	 * @param gene1 a list of exons frome gene1
	 * @param gene2 a list of exons from gene2
	 * @return return the exons that overlaps between genes ignoring strand
	 */
	public List<Feature> findOverlappingExons(List<Feature> gene1, List<Feature> gene2) {
			assert(GTFFeatureUtil.isSorted(gene1));
			assert(GTFFeatureUtil.isSorted(gene2));
			
			
			LinkedList<Feature> overlap = new LinkedList<>();
			
			if(gene1.size() == 0 || gene2.size() == 0) {
				return overlap;
			} else if(gene1.get(0).seqname().equals(gene2.get(0).seqname())) {
				overlap.addAll(findOverlappingExonsKnownOverlap(gene1,gene2));
			}
				
	
	
			return overlap;
			
	}
	
	
	public Set<Feature> findOverlappingExonsKnownOverlap(List<Feature> gene1, List<Feature> gene2) {
		assert(GTFFeatureUtil.isSorted(gene1));
		assert(GTFFeatureUtil.isSorted(gene2));
		
		Set<Feature> overlappingExons = new HashSet<>();
		
		for(Feature e1 : gene1) {
			for(Feature e2 : gene2) {
				if(e1.location().plus().overlaps(e2.location().plus())) {
					overlappingExons.add(e1);
					overlappingExons.add(e2);
				} else if(e2.location().plus().isAfter(e1.location().plus())) {
					break;
				}
			}
		}
		
		return overlappingExons;
	}
	

	
	
	
}
