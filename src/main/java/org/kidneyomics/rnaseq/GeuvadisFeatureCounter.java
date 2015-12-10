package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.List;import java.util.Map;
import java.util.Set;
import java.util.LinkedList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.FeatureComparator;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.FeatureMerger;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.SAMRecordToFeatureConverter;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import htsjdk.samtools.SAMRecord;

@Component
public class GeuvadisFeatureCounter implements FeatureCounter {

	//Chr => sorted Feature array for binary search
	private HashMap<String,Feature[]> chromosomeFeatures = new HashMap<>();
	
	HashMap<String,Feature[]> getChromosomeFeatures() {
		return this.chromosomeFeatures;
	}
	
	private HashMap<String,Integer> chromosomeLongestFeature = new HashMap<>();
	
	HashMap<String,Integer> getChromosomeLongestFeature() {
		return this.chromosomeLongestFeature;
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
	
	@Autowired
	public GeuvadisFeatureCounter(FindOverlappingFeatures findOverlappingFeatures, LoggerService loggerService) {
		this.findOverlappingFeatures = findOverlappingFeatures;
		this.logger = loggerService.getLogger(this);
	}
	
	@Override
	public void buildFeatures(File gtf, String type) throws IOException {
		
		// gene -> LinkedList<Feature>
		HashMap<String, LinkedList<Feature>> geneFeatures = new HashMap<>();
		
		if(!type.equals("exon")) {
			throw new UnsupportedOperationException("only exon counting is supported");
		}
		
		
		logger.info("Reading gtf...");
		//Read GTF
		//Store all exons for each gene
		int totalFeatures = 0;
		try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
			for(Feature feature : reader) {
				
				if(feature.type().equals("exon")) {
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
		
		
		
		HashMap<String, LinkedList<Feature>> chromosomeFeaturesUnsorted = new HashMap<>();
		
		logger.info("Merging features across each gene and inserting into separate chromosome...");
		//Merge features per gene and store
		for(Map.Entry<String, LinkedList<Feature>> entry : geneFeatures.entrySet()) {
			List<Feature> mergedFeatures = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(entry.getValue());
			
			/*
			 * create feature counts and set ids on features
			 */
			for(Feature feature : mergedFeatures) {
				FeatureCount fc = new FeatureCount(feature);
				String id = fc.getId();
				feature.getAttributes().put("id", id);
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
		FeatureComparator comparator = new FeatureComparator();
		for(Map.Entry<String, LinkedList<Feature>> entry : chromosomeFeaturesUnsorted.entrySet()) {
			Collections.sort(entry.getValue(),comparator);
		}
		
		logger.info("Converting chomosome features to array for binary search...");
		//Store as array per chromosome for binary search
		for(Map.Entry<String, LinkedList<Feature>> entry : chromosomeFeaturesUnsorted.entrySet()) {
			Feature[] featureArray = new Feature[entry.getValue().size()];
			chromosomeFeatures.put(entry.getKey(), entry.getValue().toArray(featureArray));
		}
		
		logger.info("Finding longest feature per chromosome...");
		//Find longest feature per chromosome
		for(Map.Entry<String, Feature[]> entry : chromosomeFeatures.entrySet()) {
			int longest = findOverlappingFeatures.findLongest(entry.getValue());
			chromosomeLongestFeature.put(entry.getKey(), longest);
		}
		
		
	}

	@Override
	public Collection<String> getFeatureIds() {
		return featureCounts.keySet();
	}

	@Override
	public FeatureCount getCounts(String featureId) {
		return featureCounts.get(featureId);
	}

	Set<Feature> getMappedRegionsForMate(SAMRecord mate, int longest, Feature[] features) {
		List<Feature> featuresForRecord = converter.convert(mate);
		Set<Feature> featuresForMate = new HashSet<>();
		//boolean isAmbiguous = false;
		for(Feature target : featuresForRecord) {
			//if length is greater than 1 then we are mapping to more than one gene and the read is therefore ambiguous
			//the reason for this is that the exons have been merged on a per gene level
			List<Feature> mappedFeatures = findOverlappingFeatures.findOverlappingFeatures(longest, features, target);
			featuresForMate.addAll(mappedFeatures);
			/*
			if(mappedFeatures.size() > 1) {
				isAmbiguous = true;
				break;
			} else if(mappedFeatures.size() == 1) {
				//exactly one mapped feature
				featuresForMate.add(mappedFeatures.get(0));
			}
			*/
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
	
	private void addToFeatures(double fraction, Set<Feature> features) {
		for(Feature feature : features) {
			FeatureCount fc = featureCounts.get(feature.getAttribute("id"));
			fc.addToCount(fraction);
		}
	}
	
	@Override
	public void count(SAMRecordPair samRecordPair) {
		String chr = samRecordPair.getMate1().getReferenceName();
		
		int longest = chromosomeLongestFeature.get(chr);
		Feature[] features = chromosomeFeatures.get(chr);

		/*
		 * For stranded data will have to get features for strand separately
		 */
		
		if(samRecordPair.bothPairsAligned()) {
		
			totalCount += 2;
			
			Set<Feature> mappedFeaturesAcrossChunksForMate1 = getMappedRegionsForMate(samRecordPair.getMate1(), longest, features);
			boolean isAmbiguousMate1 = mapToMultipleGenes(mappedFeaturesAcrossChunksForMate1);
			boolean isUnmappedMate1 = mappedFeaturesAcrossChunksForMate1.size() == 0;
			
			
			Set<Feature> mappedFeaturesAcrossChunksForMate2 = getMappedRegionsForMate(samRecordPair.getMate2(), longest, features);
			boolean isAmbiguousMate2 = mapToMultipleGenes(mappedFeaturesAcrossChunksForMate2);
			boolean isUnmappedMate2 = mappedFeaturesAcrossChunksForMate2.size() == 0;
			
			
			if(isAmbiguousMate1 && isAmbiguousMate2) {
				ambiguousCount +=2;
			} else if(isAmbiguousMate1) {
				ambiguousCount++;
				
				if(!isUnmappedMate2) {
					//mate2 is mapped to at least one feature
					double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate2.size();
					addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate2);
					mappedCount++;
				} else {
					unmappedCount++;
				}
				
				
			} else if(isAmbiguousMate2) {
				
				ambiguousCount++;
				
				if(!isUnmappedMate1) {
					//mate 1 is mapped to at least one feature
					double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate1.size();
					addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate1);
					mappedCount++;
				} else {
					unmappedCount++;
				}
				
			} else {
				//Both are not ambiguous individually
				
				HashSet<Feature> union = new HashSet<>();
				union.addAll(mappedFeaturesAcrossChunksForMate1);
				union.addAll(mappedFeaturesAcrossChunksForMate2);
				
				if(mapToMultipleGenes(union)) {
					ambiguousCount += 2;
				} else if(isUnmappedMate1 && isUnmappedMate2) {
					unmappedCount +=2;
				} else if(isUnmappedMate1) {
					//Mate 2 is mapped otherwise the above case
					double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate2.size();
					addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate2);
					mappedCount++;
				} else if(isUnmappedMate2) {
					//Mate 1 is mapped otherwise the above case
					double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate1.size();
					addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate1);
					mappedCount++;
				} else {
					{
						double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate1.size();
						addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate1);
					}
					{
						double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate2.size();
						addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate2);
					}
					mappedCount+=2;
				}
				
			}
			
		} else {
			//Only use first mate
			totalCount++;
			Set<Feature> mappedFeaturesAcrossChunksForMate1 = getMappedRegionsForMate(samRecordPair.getMate1(), longest, features);
			boolean isAmbiguousMate1 = mapToMultipleGenes(mappedFeaturesAcrossChunksForMate1);
			boolean isUnmappedMate1 = mappedFeaturesAcrossChunksForMate1.size() == 0;
			
			if(isAmbiguousMate1) {
				ambiguousCount++;
			} else if(isUnmappedMate1) {
				unmappedCount++;
			} else {
				double fracToAdd = 1.0 / mappedFeaturesAcrossChunksForMate1.size();
				addToFeatures(fracToAdd,mappedFeaturesAcrossChunksForMate1);
				mappedCount++;
			}
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
	
	

}
