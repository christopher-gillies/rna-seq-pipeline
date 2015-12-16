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

import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.ExonFilter;
import org.kidneyomics.gtf.FeatureComparator;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.FeatureMerger;
import org.kidneyomics.gtf.FeaturePair;
import org.kidneyomics.gtf.FindOverlappingExonsBetweenGenes;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.FindOverlappingGenePairs;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFFeatureUtil;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GeneOrExonFilter;
import org.kidneyomics.gtf.NoOpReadLogger;
import org.kidneyomics.gtf.ReadLogger;
import org.kidneyomics.gtf.RemoveRetainedIntronFilter;
import org.kidneyomics.gtf.SAMRecordToFeatureConverter;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

import htsjdk.samtools.SAMRecord;

@Scope("prototype")
@Component("gTExFeatureCounter")
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
	
	//Read logger
	private ReadLogger readLogger = new NoOpReadLogger();
	
	@Override
	public void setReadLogger(ReadLogger readLogger) {
		this.readLogger = readLogger;
	}
	
	
	//Chr => sorted Feature array for binary search
	private HashMap<String,Feature[]> chromosomeFeatures = new HashMap<>();
	
	
	HashMap<String,Feature[]> getChromosomeFeatures() {
		return this.chromosomeFeatures;
	}
	
	
	//Feature id -> FeatureCount
	private HashMap<String,FeatureCount> featureCounts = new HashMap<>();
	
	private SAMRecordToFeatureConverter converter = new SAMRecordToFeatureConverter();
	
	private double unmappedCount = 0.0;
	private double ambiguousCount = 0.0;
	private double mappedCount = 0.0;
	private double totalCount = 0.0;
	private long partiallyUnmappedReads = 0;
	
	private FindOverlappingFeatures findOverlappingFeatures;
	
	private Logger logger;
	
	private HashSet<Feature> featuresRemovedForAnalysis = null;
	
	private HashMap<String, Feature> geneInfo = new HashMap<>();
	
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

		if(!type.equals("exon")) {
			throw new UnsupportedOperationException("only exon counting is supported");
		}
		
		
		logger.info("Reading annotation gtf...");
		//Read GTF
		//Store all exons for each gene
		int totalFeatures = 0;
		int totalGenes = 0;
		int linesRead = 0;
		try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
			reader.addFilter(new GeneOrExonFilter()).addFilter(new RemoveRetainedIntronFilter());
			
						
			for(Feature feature : reader) {
				linesRead++;
				if(linesRead % 100000 == 0) {
					logger.info("Lines read: " + linesRead);
					logger.info("Last read feature: " + GTFFeatureRenderer.render(feature));
				}
				if(feature.type().equals("gene")) {
					geneInfo.put(feature.getAttribute("gene_id"), feature);
					totalGenes++;
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
		
		logger.info("Total genes " + totalGenes);
		logger.info("Total exons " + totalFeatures);
		
		/*
		 * Checking GTF Validity
		 */
		logger.info("Validating that every exon has at least some gene entry from GTF file");
		for(String key : geneFeatures.keySet()) {
			if(!geneInfo.containsKey(key)) {
				throw new IllegalStateException(key + " does not have a gene entry in the GTF");
			}
		}
		logger.info("All gene ids found for exons");
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
		
		logger.info("Finding overlapping exons between genes");
		HashSet<Feature> geneOverlapToRemove = new HashSet<>();
		FindOverlappingGenePairs findOverlappingGenePairs = new FindOverlappingGenePairs(geneInfo.values());
		List<FeaturePair> overlappingGenes = findOverlappingGenePairs.getOverlappingGenes();
		for(FeaturePair pair : overlappingGenes) {
			//logger.info("First of pair: " + pair.getFirst().getAttribute("gene_id"));
			//logger.info("Second of pair: " + pair.getSecond().getAttribute("gene_id"));
			List<Feature> gene1 = geneFeatures.get(pair.getFirst().getAttribute("gene_id"));
			List<Feature> gene2 = geneFeatures.get(pair.getSecond().getAttribute("gene_id"));
			//Some genes may have only retained introns...
			if(gene1 != null && gene1.size() > 0 && gene2 != null && gene2.size() > 0) {
				List<Feature> overlappingExons = findOverlappingExonsBetweenGenes.findOverlappingExons(gene1,gene2);
				if(overlappingExons.size() > 0) {
					geneOverlapToRemove.addAll(overlappingExons);
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
			//logger.info("GENE " + geneKey + " has " + features.size());
			Iterator<Feature> iter = features.iterator();
			while(iter.hasNext()) {
				Feature f = iter.next();
				if(geneOverlapToRemove.contains(f)) {
					iter.remove();
					countRemovedFeatures++;
				}
			}
			//logger.info("GENE " + geneKey + " has " + features.size() + " after removing overlapping exons");
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
	
	static List<Feature> computeGeneLevelExpression(Map<String,FeatureCount> featureCounts, Map<String,Feature> geneInfos, double numberOfReads) {
		List<Feature> geneFeatures = new LinkedList<Feature>();
		
		HashMap<String,List<FeatureCount>> geneFeaturesMap = new HashMap<>();
		
		//Organize exons by gene
		Collection<FeatureCount> exons = featureCounts.values();
		for(FeatureCount fc : exons) {
			String geneId = fc.getFeature().getAttribute("gene_id");
			if(geneFeaturesMap.containsKey(geneId)) {
				geneFeaturesMap.get(geneId).add(fc);
			} else {
				LinkedList<FeatureCount> fcs = new LinkedList<>();
				fcs.add(fc);
				geneFeaturesMap.put(geneId, fcs);
			}
		}
		
		//compute total counts and total length
		for(Map.Entry<String, List<FeatureCount>> entry : geneFeaturesMap.entrySet()) {
			List<FeatureCount> fcsForGene = entry.getValue();
			Collections.sort(fcsForGene);
			Feature geneInfo = geneInfos.get(entry.getKey());
			
			
			if(geneInfo == null) {
				throw new IllegalStateException(entry.getKey() + " has no entry in geneInfos");
			}
			
			if(!geneInfo.type().equals("gene")) {
				throw new IllegalStateException(geneInfo + " is not a gene");
			}
			
			//get total length and counts
			double count = 0.0;
			int length = 0;
			List<String> exonIds = new LinkedList<>();
			for(FeatureCount fc : fcsForGene) {
				count += fc.getCount();
				length += fc.getFeature().location().length();
				exonIds.add(fc.getId());
			}
			
			//comput rpkm
			double rpkm = count / length / numberOfReads * Math.pow(10, 9);
			
			HashMap<String,String> info = new HashMap<>();
			info.put("reads", Double.toString(count));
			info.put("RPKM", Double.toString(rpkm));
			info.put("length", Integer.toString(length));
			info.put("exons", StringUtils.join(exonIds, ","));
			
			Feature newGeneInfo = GTFFeatureBuilder.addAttributesToFeature(geneInfo, info);
			
			geneFeatures.add(newGeneInfo);
		}
		return geneFeatures;
	}
	

	static Map<Feature,Integer> getMappedRegionsForMate(SAMRecord mate, Feature[] features, SAMRecordToFeatureConverter converter, FindOverlappingFeatures findOverlappingFeatures) {
		List<Feature> featuresForRecord = converter.convert(mate);
		HashMap<Feature,Integer> featuresForMate = new HashMap<>();
		//boolean isAmbiguous = false;
		for(Feature target : featuresForRecord) {
			//logger.info(GTFFeatureRenderer.render(target));
			//we can assume the result is non overlapping
			List<Feature> mappedFeatures = findOverlappingFeatures.findOverlappingFeatures( features, target);
			//logger.info("Number of mapped features: " + mappedFeatures.size());
			assert(GTFFeatureUtil.hasNoOverlapIgnoreStrand(mappedFeatures));
			for(Feature mapped : mappedFeatures) {
				//logger.info(GTFFeatureRenderer.render(mapped));
				if(target.location().plus().overlaps(mapped.location().plus())) {
					int overlap = target.location().plus().intersection(mapped.location().plus()).length();
					//logger.info("overlap: " + overlap);
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
	
	
	static boolean mapToMultipleGenes(Set<Feature> features) {
		HashSet<String> geneIds = new HashSet<>();
		for(Feature feature : features) {
			String geneId = feature.getAttribute("gene_id");
			geneIds.add(geneId);
		}
		return geneIds.size() > 1;
	}
	
	static double addToFeatures(int totalMappedBases, Map<Feature,Integer> features, Map<String,FeatureCount> featureCounts) {
		assert(features != null);
		assert(featureCounts != null);
		
		if(features.size() == 0) {
			return 1.0;
		}
		
		/*
		 * Count the total mapping
		 * just going to use total mapped bases
		 */
		//double total = 0.0;
		//for(Map.Entry<Feature, Integer> entry : features.entrySet()) {
		//	total += entry.getValue();
		//}
		double mapped = 0.0;
		for(Map.Entry<Feature, Integer> entry : features.entrySet()) {
			Feature feature = entry.getKey();
			double count = entry.getValue();
			mapped += count;
			String id = feature.getAttribute("id");
			if(id == null) {
				throw new IllegalArgumentException("Feature must have id attribute\n" + GTFFeatureRenderer.render(feature));
			}
			FeatureCount fc = featureCounts.get(id);
			if(fc == null) {
				throw new IllegalStateException("No FeatureCount object found for id: " + id);
			}
			
			fc.addToCount(count / (double) totalMappedBases);
		}
		double result = ((double) totalMappedBases - mapped)   /  (double) totalMappedBases;
		if(result < -0.0001) {
			throw new IllegalStateException("There was an issue counting reads...\nresult=" + result + "\ntotalMappedBases=" + totalMappedBases + "\nmapped=" + mapped);
		}
		return result;
	}
	
	static Map<Feature,Integer> unionFeaturesForMates(Map<Feature,Integer> mappedFeaturesAcrossChunksForMate1, Map<Feature,Integer> mappedFeaturesAcrossChunksForMate2) {
		//This map will give us the number of bases overlapping each interval, thus we can fractionally count
		Map<Feature,Integer> union = new HashMap<>();
		
		
		for(Map.Entry<Feature, Integer> entry : mappedFeaturesAcrossChunksForMate1.entrySet()) {
			if(union.containsKey(entry.getKey())) {
				int currentVal = union.get(entry.getKey());
				union.put(entry.getKey(), entry.getValue() + currentVal);
			} else {
				union.put(entry.getKey(), entry.getValue());
			}
		}
		
		for(Map.Entry<Feature, Integer> entry : mappedFeaturesAcrossChunksForMate2.entrySet()) {
			if(union.containsKey(entry.getKey())) {
				int currentVal = union.get(entry.getKey());
				union.put(entry.getKey(), entry.getValue() + currentVal);
			} else {
				union.put(entry.getKey(), entry.getValue());
			}
		}
		
		return union;
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
			
			Map<Feature,Integer> mappedFeaturesAcrossChunksForMate1 = getMappedRegionsForMate(samRecordPair.getMate1(), features, converter, findOverlappingFeatures);	
			boolean isUnmappedMate1 = mappedFeaturesAcrossChunksForMate1.size() == 0;
			int mate1MappedBases = converter.getNumberOfMappedBases(samRecordPair.getMate1());
			
			Map<Feature,Integer> mappedFeaturesAcrossChunksForMate2 = getMappedRegionsForMate(samRecordPair.getMate2(), features, converter, findOverlappingFeatures);
			boolean isUnmappedMate2 = mappedFeaturesAcrossChunksForMate2.size() == 0;
			int mate2MappedBases = converter.getNumberOfMappedBases(samRecordPair.getMate2());
			
			Map<Feature,Integer> union = unionFeaturesForMates(mappedFeaturesAcrossChunksForMate1,mappedFeaturesAcrossChunksForMate2);
			
			int totalMappedBases = mate1MappedBases + mate2MappedBases;
			
			
			
			//logger.info("Ambiguous count: " + ambiguousCount);
			//logger.info("Read name: " + samRecordPair.getMate1().getReadName());
			//logger.info("Coordinates: " + samRecordPair.getMate1().getReferenceName() + ":" + samRecordPair.getMate1().getAlignmentStart() + "-" + samRecordPair.getMate1().getAlignmentEnd());
			if(mapToMultipleGenes(union.keySet())) {
				ambiguousCount += 1;
				readLogger.logRead("ambiguous", samRecordPair.getMate1(),mappedFeaturesAcrossChunksForMate1);
				readLogger.logRead("ambiguous", samRecordPair.getMate2(),mappedFeaturesAcrossChunksForMate2);
			} else if(isUnmappedMate1 && isUnmappedMate2) {
				unmappedCount += 1;
				readLogger.logRead("unmapped", samRecordPair.getMate1(),mappedFeaturesAcrossChunksForMate1);
				readLogger.logRead("unmapped", samRecordPair.getMate2(),mappedFeaturesAcrossChunksForMate2);
			} else if(isUnmappedMate1) {
				//Mate 2 is mapped otherwise the above case
				double unmappedFrac = addToFeatures(totalMappedBases, mappedFeaturesAcrossChunksForMate2, featureCounts);
				mappedCount += (1 - unmappedFrac);
				unmappedCount += unmappedFrac;
				if(Math.abs(unmappedFrac) >= 0.0001) {
					partiallyUnmappedReads++;
					readLogger.logRead("partially_unmapped", samRecordPair.getMate1(),mappedFeaturesAcrossChunksForMate1);
					readLogger.logRead("partially_unmapped", samRecordPair.getMate2(),mappedFeaturesAcrossChunksForMate2);
				} 
			} else if(isUnmappedMate2) {
				//Mate 1 is mapped otherwise the above case
				double unmappedFrac = addToFeatures(totalMappedBases, mappedFeaturesAcrossChunksForMate1, featureCounts);
				mappedCount += (1 - unmappedFrac);
				unmappedCount += unmappedFrac;
				if(Math.abs(unmappedFrac) >= 0.0001) {
					partiallyUnmappedReads++;
					readLogger.logRead("partially_unmapped", samRecordPair.getMate1(),mappedFeaturesAcrossChunksForMate1);
					readLogger.logRead("partially_unmapped", samRecordPair.getMate2(),mappedFeaturesAcrossChunksForMate2);
				} 
			} else {
				double unmappedFrac = addToFeatures(totalMappedBases, union, featureCounts);
				mappedCount += (1 - unmappedFrac);
				unmappedCount += unmappedFrac;
				if(Math.abs(unmappedFrac) >= 0.0001) {
					partiallyUnmappedReads++;
					readLogger.logRead("partially_unmapped", samRecordPair.getMate1(),mappedFeaturesAcrossChunksForMate1);
					readLogger.logRead("partially_unmapped", samRecordPair.getMate2(),mappedFeaturesAcrossChunksForMate2);
				} 
			}
			//logger.info("Ambiguous count: " + ambiguousCount);	
			
			
		} else {
			throw new IllegalStateException("Read is not properly paired. Because its pair was not found. " + samRecordPair.getMate1().getReadName());
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
		
		if(!this.validState()) {
			logInfo();
			throw new IllegalStateException("Counts do not add up correctly. Last read: " + samRecordPair.getMate1().getReadName() + " " + samRecordPair.getMate1().getReferenceName() + ":" + samRecordPair.getMate1().getAlignmentStart());
		}
		
	}
	
	@Override
	public double getAmbiguousReadCount() {
		return this.ambiguousCount;
	}

	@Override
	public double getUnmappedReadCount() {
		return this.unmappedCount;
	}

	@Override
	public double getMappedReadCount() {
		return this.mappedCount;
	}

	@Override
	public double getTotalCount() {
		return this.totalCount;
	}

	@Override
	public List<Feature> getCounts() {
		List<Feature> newFeatures = new LinkedList<>();
		for(FeatureCount fc : this.featureCounts.values()) {
			Feature f = fc.getFeature();
			
			Map<String,String> newAtts = new HashMap<>();
			newAtts.put("id", fc.getId());
			newAtts.put("reads", Double.toString(fc.getCount()));
			newAtts.put("RPKM", Double.toString(fc.getRPKM(this.getTotalCount())));
			
			Feature newF = GTFFeatureBuilder.addAttributesToFeature(f, newAtts);
			newFeatures.add(newF);
		}
		return newFeatures;
	}

	@Override
	public long getNumberOfPartiallyUnmappedReads() {
		return this.partiallyUnmappedReads;
	}

	@Override
	public boolean validState() {
		return Math.abs( this.totalCount - (this.mappedCount + this.unmappedCount + this.ambiguousCount)) < 0.01; 
	}

	@Override
	public void logInfo() {
		logger.info("Current status:");
		logger.info("Total counts reads: " + this.getTotalCount());
		logger.info("Total mapped read count: " + this.getMappedReadCount());
		logger.info("Total unmapped read count: " + this.getUnmappedReadCount());
		logger.info("Total ambiguous read count: " + this.getAmbiguousReadCount());
		logger.info("Partially unmapped reads: " + this.getNumberOfPartiallyUnmappedReads());
	}

	@Override
	public ReadLogger getReadLogger() {
		return this.readLogger;
	}

	@Override
	public List<Feature> getGeneCounts() {
		return GTExFeatureCounter.computeGeneLevelExpression(this.featureCounts, this.geneInfo, this.getTotalCount());
	}
	
	

}
