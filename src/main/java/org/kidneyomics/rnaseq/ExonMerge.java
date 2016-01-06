package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.GTFFeatureRenderer;
import org.kidneyomics.gtf.GTFReader;
import org.kidneyomics.gtf.GeneOrExonFilter;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Component;

@Component
@Scope("prototype")
public class ExonMerge implements ApplicationCommand {

	
	private ApplicationOptions applicationOptions;
	
	private Logger logger;
	
	private SampleGTFReader sampleGtfReader;
	
	private QuantificationFactory quantificationFactory;
	
	/*
	 * Create datastores
	 */
	HashMap<String,MutableQuantification> exonCounts = new HashMap<>();
	HashMap<String,MutableQuantification> exonRpkm = new HashMap<>();
	HashMap<String,MutableQuantification> geneCounts = new HashMap<>();
	HashMap<String,MutableQuantification> geneRpkm = new HashMap<>();
	
	
	@Autowired
	ExonMerge(LoggerService loggerService, ApplicationOptions applicationOptions, SampleGTFReader sampleGtfReader, QuantificationFactory quantificationFactory) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
		this.sampleGtfReader = sampleGtfReader;
		this.quantificationFactory = quantificationFactory;
	}
	
	public void writeOutMatrices() throws IOException {
		String gtfList = applicationOptions.getFileIn();
		String outDir = applicationOptions.getOutputDirectory();
		
		
		/*
		 * Read sample information
		 */
		logger.info("Reading sample gtf file list");
		List<SampleGTF> sampleGtfs = sampleGtfReader.readFile(new File(gtfList));
		
		if(sampleGtfs.size() == 0) {
			throw new IllegalArgumentException("Please specify make sure the list of sample GTF files is not empty");
		}
		
		logger.info("Sorting samples by id");
		Collections.sort(sampleGtfs);
		List<String> sampleIds = sampleGtfReader.getIds(sampleGtfs);
		Collections.sort(sampleIds);
		
		assert(sampleGtfs.size() > 0);
		Iterator<SampleGTF> iter = sampleGtfs.iterator();
		
		//loop through all the samples and read their gtfs and store counts
		logger.info("Reading sample gtfs");
		while(iter.hasNext()) {
			//get next sample
			SampleGTF sgtf = iter.next();
			String sampleId = sgtf.getSampleId();
			
			//get file
			File gtf = sgtf.getFile();
			
			logger.info("Reading sample: " + sampleId + " " + gtf.getAbsolutePath());
			
			//Create reader
			try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
				reader.addFilter(new GeneOrExonFilter());
				
				//read through lines of gtf
				for(Feature feature : reader) {
					
					//organize by type
					if(feature.type().equals("exon")) {
						String featureId = feature.getAttribute("id");
						if(featureId == null) {
							throw new IllegalStateException("This feature has no exon id: " + GTFFeatureRenderer.render(feature));
						}
						//count exon reads
						countFeature(sampleId, sampleIds, featureId,feature,exonCounts,QuantificationType.EXON_COUNTS);
						//count exon rpkm
						countFeature(sampleId, sampleIds, featureId,feature,exonRpkm,QuantificationType.EXON_RPKM);
					} else if(feature.type().equals("gene")) {
						String featureId = feature.getAttribute("gene_id");
						//count gene reads
						countFeature(sampleId, sampleIds, featureId,feature,geneCounts,QuantificationType.GENE_COUNTS);
						//count gene rpkm
						countFeature(sampleId, sampleIds, featureId,feature,geneRpkm,QuantificationType.GENE_RPKM);
					} else {
						throw new IllegalStateException("Error read a nonexon or nongene entry for " + gtf.getAbsolutePath());
					}
				}
			}
		}
		
		logger.info("Read all samples");
		
		//write results
		logger.info("Writing out results");
		
		writeResults(exonCounts,outDir + "/exon.counts.txt");
		writeResults(exonRpkm,outDir + "/exon.rpkm.txt");
		writeResults(geneCounts,outDir + "/gene.counts.txt");
		writeResults(geneRpkm,outDir + "/gene.rpkm.txt");
		
		logger.info("Complete");
	}
	
	private void writeResults(HashMap<String,MutableQuantification> map, String outfile) throws IOException {
		ArrayList<MutableQuantification> mqs = new ArrayList<>(map.values().size());
		mqs.addAll(map.values());
		Collections.sort(mqs);
		QuantificationUtils.writeQuantificationMatrix(mqs, outfile);
	}
	
	private void countFeature(String sampleId, List<String> sampleIds, String featureId, Feature feature, HashMap<String,MutableQuantification> map, QuantificationType type) {
		MutableQuantification mq = null;
		//Check if the feature exists
		//if not then create it
		if(!map.containsKey(featureId)) {
			mq = quantificationFactory.getQuantification(feature, sampleIds, type);
			map.put(featureId, mq);
		} else {
			mq = map.get(featureId);
		}
		//count the feature
		switch(type) {
		case EXON_COUNTS:
		case GENE_COUNTS:
			mq.putSampleCount(sampleId, feature);
			break;
		case EXON_RPKM:
		case GENE_RPKM:
			mq.putSampleRPKM(sampleId, feature);
			break;
			default:
				throw new IllegalStateException(type + " not supported");
		}
	}

	@Override
	public void doWork() throws Exception {
		writeOutMatrices();
	}
	
	
}
