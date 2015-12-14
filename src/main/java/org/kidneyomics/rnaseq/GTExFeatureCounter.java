package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.List;import java.util.Map;
import java.util.Set;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.ExonFilter;
import org.kidneyomics.gtf.FeatureComparator;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.FeatureMerger;
import org.kidneyomics.gtf.FindOverlappingExonsBetweenGenes;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFFeatureUtil;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GeneOrExonFilter;
import org.kidneyomics.gtf.RemoveRetainedIntronFilter;
import org.kidneyomics.gtf.SAMRecordToFeatureConverter;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import htsjdk.samtools.SAMRecord;

@Component
public class GTExFeatureCounter implements FeatureCounter {

	/*
	 * 
	 * Note: count read pairs
	 * 
	 * Unqiuely mappped reads
	 * properly paired
	 * maybe we will do the edit distance thing
	 * 
	 * (1) merge overlapping exons per gene
	 * (2) remove retained introns
	 * (3) Discard intervals corresponding to multiple genes
	 * 
	 */
	
	//Chr => sorted Feature array for binary search
	private HashMap<String,Feature[]> chromosomeFeatures = new HashMap<>();
	
	
	HashMap<String,Feature[]> getChromosomeFeatures() {
		return this.chromosomeFeatures;
	}
	
	
	//Feature id -> FeatureCount
	private HashMap<String,FeatureCount> featureCounts = new HashMap<>();
	
	private SAMRecordToFeatureConverter converter = new SAMRecordToFeatureConverter();
	
	private long unmappedCount = 0;
	private long ambiguousCount = 0;
	private long mappedCount = 0;
	private long totalCount = 0;
	
	
	private FindOverlappingFeatures findOverlappingFeatures;
	
	private Logger logger;
	
	private HashSet<Feature> featuresRemovedForAnalysis = null;
	
	public HashSet<Feature> getRemovedFeatures() {
		return this.featuresRemovedForAnalysis;
	}
	
	@Autowired
	public GTExFeatureCounter(FindOverlappingFeatures findOverlappingFeatures, LoggerService loggerService) {
		this.findOverlappingFeatures = findOverlappingFeatures;
		this.logger = loggerService.getLogger(this);
	}
	
	@Override
	public void buildFeatures(File gtf, String type) throws IOException {
		
		// gene -> LinkedList<Feature>
		HashMap<String, LinkedList<Feature>> geneFeatures = new HashMap<>();
		HashMap<String, Feature> geneInfo = new HashMap<>();
		
		if(!type.equals("exon")) {
			throw new UnsupportedOperationException("only exon counting is supported");
		}
		
		
		logger.info("Reading gtf...");
		//Read GTF
		//Store all exons for each gene
		int totalFeatures = 0;
		try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
			reader.addFilter(new GeneOrExonFilter()).addFilter(new RemoveRetainedIntronFilter());
			
						
			for(Feature feature : reader) {
				if(feature.type().equals("gene")) {
					geneInfo.put(feature.getAttribute("gene_id"), feature);
				} else if(feature.type().equals("exon")) {
					totalFeatures++;
					String geneId = feature.getAttribute("gene_id");
					if(geneFeatures.containsKey(geneId)) {
						geneFeatures.get(geneId).add(feature);
					} else {
						LinkedList<Feature> features = new LinkedList<Feature>();
						features.add(feature);
						geneFeatures.put(geneId, features);
					}
				}
				
				
				
			}
		}
		
		logger.info("Total exons " + totalFeatures);
		

		/*
		 * sorting exons per gene
		 */
		
		logger.info("Sorting features per gene...");
		//Sort features per chromosome
		for(Map.Entry<String, LinkedList<Feature>> entry : geneFeatures.entrySet()) {
			GTFFeatureUtil.sortFeatures(entry.getValue());
		}
		
		/*
		 * find overlapping features between genes
		 */
		ArrayList<String> geneKeys = new ArrayList<String>(geneFeatures.keySet().size());
		geneKeys.addAll(geneFeatures.keySet());
		
		FindOverlappingExonsBetweenGenes findOverlappingExonsBetweenGenes = new FindOverlappingExonsBetweenGenes();
		
		/***
		 * THIS CODE BELOW MAY BE A BIG BOTTLENECK BUT WE WILL SEE
		 */
		logger.info("Finding overlapping exons between genes");
		HashSet<Feature> geneOverlapToRemove = new HashSet<>();
		for(int i = 0; i < geneKeys.size() - 1; i++) {
			for(int j = i + 1; j < geneKeys.size(); j++) {
				String key1 = geneKeys.get(i);
				String key2 = geneKeys.get(j);
				
				if(GTFFeatureUtil.overlapIgnoreStrand(geneInfo.get(key1), geneInfo.get(key2))) {
					List<Feature> gene1 = geneFeatures.get(key1);
					List<Feature> gene2 = geneFeatures.get(key2);
					List<Feature> overlappingExons = findOverlappingExonsBetweenGenes.findOverlappingExons(gene1,gene2);
					if(overlappingExons.size() > 0) {
						geneOverlapToRemove.addAll(overlappingExons);
					}
				}
			}
		}
		
		logger.info("Storing removed exons");
		featuresRemovedForAnalysis = geneOverlapToRemove;
		logger.info("Total overlapping features between genes: " + geneOverlapToRemove.size());
		
		/*
		 * remove overlapping features between genes
		 */
		int countRemovedFeatures = 0;
		logger.info("Removing overlapping exons per gene");
		for(String geneKey : geneKeys) {
			List<Feature> features = geneFeatures.get(geneKey);
			logger.info("GENE " + geneKey + " has " + features.size());
			Iterator<Feature> iter = features.iterator();
			while(iter.hasNext()) {
				Feature f = iter.next();
				if(geneOverlapToRemove.contains(f)) {
					iter.remove();
					countRemovedFeatures++;
				}
			}
			logger.info("GENE " + geneKey + " has " + features.size() + " after removing overlapping exons");
		}
		logger.info("Total overlapping features actually removed: " + countRemovedFeatures);
		
		HashMap<String, LinkedList<Feature>> chromosomeFeaturesUnsorted = new HashMap<>();
		
		logger.info("Merging features across each gene and inserting into separate chromosome...");
		//Merge features per gene and store
		for(Map.Entry<String, LinkedList<Feature>> entry : geneFeatures.entrySet()) {
			
			if(entry.getValue().size() == 0) {
				continue;
			}
			
			//GTFFeatureUtil.sortFeatures(entry.getValue());
			List<Feature> mergedFeatures = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(entry.getValue());
			
			/*
			 * create feature counts and set ids on features
			 */
			for(Feature feature : mergedFeatures) {
				FeatureCount fc = new FeatureCount(feature);
				String id = fc.getId();
				featureCounts.put(id, fc);
			}
			
			
			Feature first = mergedFeatures.get(0);
			String chr = first.seqname();
			/*
			 * For stranded data add + or - to end of chr
			 */
			
			if(chromosomeFeaturesUnsorted.containsKey(chr)) {
				chromosomeFeaturesUnsorted.get(chr).addAll(mergedFeatures);
			} else {
				LinkedList<Feature> chrFeatures = new LinkedList<>();
				chrFeatures.addAll(mergedFeatures);
				chromosomeFeaturesUnsorted.put(chr, chrFeatures);
			}
		}
		
		logger.info("Sorting features per chromosome...");
		//Sort features per chromosome
		for(Map.Entry<String, LinkedList<Feature>> entry : chromosomeFeaturesUnsorted.entrySet()) {
			GTFFeatureUtil.sortFeatures(entry.getValue());
		}
		
		logger.info("Validating that there are no overlapping intervals");
		
		for(Map.Entry<String, LinkedList<Feature>> entry : chromosomeFeaturesUnsorted.entrySet()) {
			logger.info("Validating that " + entry.getKey() + " has no overlap");
			if(!GTFFeatureUtil.hasNoOverlapIgnoreStrand(entry.getValue())) {
				throw new IllegalStateException("Some intervals have overlap... There was a problem in the merging process");
			}
		}
		
		logger.info("Converting chomosome features to array for binary search...");
		//Store as array per chromosome for binary search
		for(Map.Entry<String, LinkedList<Feature>> entry : chromosomeFeaturesUnsorted.entrySet()) {
			Feature[] featureArray = new Feature[entry.getValue().size()];
			chromosomeFeatures.put(entry.getKey(), entry.getValue().toArray(featureArray));
			
		}
		
		

		
	}

	@Override
	public Collection<String> getFeatureIds() {
		return Collections.unmodifiableCollection(featureCounts.keySet());
	}

	@Override
	public FeatureCount getCounts(String featureId) {
		return featureCounts.get(featureId);
	}
	
	

	Map<Feature,Integer> getMappedRegionsForMate(SAMRecord mate, Feature[] features) {
		List<Feature> featuresForRecord = converter.convert(mate);
		HashMap<Feature,Integer> featuresForMate = new HashMap<>();
		//boolean isAmbiguous = false;
		for(Feature target : featuresForRecord) {
			//we can assume the result is non overlapping
			List<Feature> mappedFeatures = findOverlappingFeatures.findOverlappingFeatures( features, target);
			assert(GTFFeatureUtil.hasNoOverlapIgnoreStrand(mappedFeatures));
			for(Feature mapped : mappedFeatures) {
				if(target.location().plus().overlaps(mapped.location().plus())) {
					int overlap = target.location().plus().intersection(mapped.location().plus()).length();
					if(featuresForMate.containsKey(mapped)) {
						int currentVal = featuresForMate.get(mapped);
						featuresForMate.put(mapped, overlap + currentVal);
					} else {
						featuresForMate.put(mapped, overlap);
					}
					
				}
			}
			
		}
		return featuresForMate;
	}
	
	
	boolean mapToMultipleGenes(Set<Feature> features) {
		HashSet<String> geneIds = new HashSet<>();
		for(Feature feature : features) {
			String geneId = feature.getAttribute("gene_id");
			geneIds.add(geneId);
		}
		return geneIds.size() > 1;
	}
	
	static void addToFeatures(Map<Feature,Integer> features, Map<String,FeatureCount> featureCounts) {
		assert(features != null);
		assert(featureCounts != null);
		
		if(features.size() == 0) {
			return;
		}
		
		/*
		 * Count the total mapping
		 */
		double total = 0.0;
		for(Map.Entry<Feature, Integer> entry : features.entrySet()) {
			total += entry.getValue();
		}
		
		for(Map.Entry<Feature, Integer> entry : features.entrySet()) {
			Feature feature = entry.getKey();
			double count = entry.getValue();
			
			String id = feature.getAttribute("id");
			if(id == null) {
				throw new IllegalArgumentException("Feature must have id attribute\n" + GTFFeatureRenderer.render(feature));
			}
			FeatureCount fc = featureCounts.get(id);
			if(fc == null) {
				throw new IllegalStateException("No FeatureCount object found for id: " + id);
			}
			
			fc.addToCount(count / total);
		}
	}
	
	@Override
	public void count(SAMRecordPair samRecordPair) {
		String chr = samRecordPair.getMate1().getReferenceName();
		
		Feature[] features = chromosomeFeatures.get(chr);

		/*
		 * For stranded data will have to get features for strand separately
		 */
		
		if(samRecordPair.bothPairsAligned()) {
		
			totalCount += 1;
			
			Map<Feature,Integer> mappedFeaturesAcrossChunksForMate1 = getMappedRegionsForMate(samRecordPair.getMate1(), features);	
			boolean isUnmappedMate1 = mappedFeaturesAcrossChunksForMate1.size() == 0;
			
			
			Map<Feature,Integer> mappedFeaturesAcrossChunksForMate2 = getMappedRegionsForMate(samRecordPair.getMate2(), features);
			boolean isUnmappedMate2 = mappedFeaturesAcrossChunksForMate2.size() == 0;
			
			Map<Feature,Integer> union = new HashMap<>();
			
			//This map will give us the number of bases overlapping each interval, thus we can fractionally count
			for(Map.Entry<Feature, Integer> entry : mappedFeaturesAcrossChunksForMate1.entrySet()) {
				if(union.containsKey(entry.getKey())) {
					int currentVal = union.get(entry.getKey());
					union.put(entry.getKey(), entry.getValue() + currentVal);
				} else {
					union.put(entry.getKey(), entry.getValue());
				}
			}
			
			
			if(mapToMultipleGenes(union.keySet())) {
				ambiguousCount += 1;
			} else if(isUnmappedMate1 && isUnmappedMate2) {
				unmappedCount += 1;
			} else if(isUnmappedMate1) {
				//Mate 2 is mapped otherwise the above case
				addToFeatures(mappedFeaturesAcrossChunksForMate2, featureCounts);
				mappedCount++;
			} else if(isUnmappedMate2) {
				//Mate 1 is mapped otherwise the above case
				addToFeatures(mappedFeaturesAcrossChunksForMate1, featureCounts);
				mappedCount++;
			} else {
				addToFeatures(union, featureCounts);
				mappedCount++;
			}
				
			
			
		} else {
			throw new IllegalStateException("Should not get here");
			//Only use first mate
//			totalCount++;
//			Map<Feature,Integer> mappedFeaturesAcrossChunksForMate1 = getMappedRegionsForMate(samRecordPair.getMate1() , features);
//			boolean isAmbiguousMate1 = mapToMultipleGenes(mappedFeaturesAcrossChunksForMate1.keySet());
//			boolean isUnmappedMate1 = mappedFeaturesAcrossChunksForMate1.size() == 0;
//			
//			if(isAmbiguousMate1) {
//				ambiguousCount++;
//			} else if(isUnmappedMate1) {
//				unmappedCount++;
//			} else {
//				double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate1.size();
//				//addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate1);
//				mappedCount++;
//			}
		}
		
	}
	
	@Override
	public long getAmbiguousReadCount() {
		return this.ambiguousCount;
	}

	@Override
	public long getUnmappedReadCount() {
		return this.unmappedCount;
	}

	@Override
	public long getMappedReadCount() {
		return this.mappedCount;
	}

	@Override
	public long getTotalCount() {
		return this.totalCount;
	}

	@Override
	public List<FeatureCount> getCounts() {
		ArrayList<FeatureCount> al = new ArrayList<FeatureCount>();
		al.addAll(this.featureCounts.values());
		return al;
	}
	
	

}
