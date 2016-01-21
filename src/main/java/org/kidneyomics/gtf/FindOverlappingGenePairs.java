package org.kidneyomics.gtf;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;

/**
 * The purpose of this program is to get a list candidate genes that overlap each other
 * @author cgillies
 *
 */
public class FindOverlappingGenePairs {
	
	private final HashMap<String,List<Feature>> genes;
	
	public FindOverlappingGenePairs(Collection<Feature> genes) {
		this.genes = new HashMap<>();
		for(Feature f : genes) {
			String chr = f.seqname();
			if(this.genes.containsKey(chr)) {
				this.genes.get(chr).add(f);
			} else {
				List<Feature> chrGenes = new ArrayList<>(10000);
				chrGenes.add(f);
				this.genes.put(chr, chrGenes);
			}
		}
		
		//sort genes
		for(Map.Entry<String, List<Feature>> entry : this.genes.entrySet()) {
			GTFFeatureUtil.sortFeatures(entry.getValue());
		}
	}
	
	
	
	public List<FeaturePair> getOverlappingGenes() {
		LinkedList<FeaturePair> pairs = new LinkedList<>();
		
		//find overlapping genes per chr
		for(Map.Entry<String, List<Feature>> entry : this.genes.entrySet()) {
			List<Feature> genesForChr = entry.getValue();
			
			//since genes are sorted we know if we find a gene that does not overlap with the current gene, gene1 then no more could overlap it
			// let G = g1, g2, g3, ... gn be the set of genes sorted by start coordinates
			// claim if gi does not overlap gj, then there is no gk such that k > j that overlaps gi. where i < j
			// suppose not, then there exists a gk where k > j such that gi does not overlap gj but gi overlaps gk
			// but if gi does not overlap gj therefore gj.start > gi.end, otherwise they would overlap. 
			// from above we have gk.start > gj.start > gi.end therefore we have a contradiction b/c gi.start < gk.start and gk.start > gi.end but there is overlap which is not possible
			for(int i = 0; i < genesForChr.size(); i++) {
				Feature gene1 = genesForChr.get(i);
				for(int j = i + 1; j < genesForChr.size(); j++) {
					Feature gene2 = genesForChr.get(j);
					if(GTFFeatureUtil.overlapIgnoreStrand(gene1, gene2)) {
						pairs.add( new FeaturePairImpl(gene1, gene2));
					} else {
						break;
					}
				}
			}
		}
		
		
		return pairs;
	}
	
	
	
	private class FeaturePairImpl implements FeaturePair {

		private final Feature first;
		private final Feature second;
		
		public FeaturePairImpl(Feature first, Feature second) {
			this.first = first;
			this.second = second;
		}
		
		@Override
		public Feature getFirst() {
			return this.first;
		}

		@Override
		public Feature getSecond() {
			return this.second;
		}

		@Override
		public boolean overlaps() {
			return this.first.seqname().equals(second.seqname()) && this.first.location().plus().overlaps(second.location().plus());
		}
		
	}
	
}
