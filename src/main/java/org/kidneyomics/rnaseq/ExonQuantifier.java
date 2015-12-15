package org.kidneyomics.rnaseq;

import java.io.File;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFWriter;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

@Component
public class ExonQuantifier {

	
	@Autowired
	ApplicationOptions applicationOptions;
	
	@Autowired
	FeatureCounter gTExFeatureCounter;
	
	Logger logger;
	
	@Autowired
	ExonQuantifier(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	
	public void quantify() throws Exception {
		
		String gtfFile = applicationOptions.getGtf();
		String fileOut = applicationOptions.getFileOut();
		String fileIn = applicationOptions.getFileIn();
		
		//Read GTF file
		//FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		//FeatureCounter featureCounter = new GTExFeatureCounter(findOverlappingFeatures, loggerService)
		
		gTExFeatureCounter.buildFeatures(new File(gtfFile), "exon");
		
		
		//Process BAM
		File fin = new File(fileIn);
		int pairsCounted = 0;
		try(BAMProcessor processor = BAMProcessor.getBAMProcessor(fin)) {
			
			SAMRecordPair pair = null;
			while( ( pair = processor.getNextReadPair()) != null) {
				gTExFeatureCounter.count(pair);
				pairsCounted++;
				if(pairsCounted % 10000 == 0) {
					logger.info("Read pairs counted: " + pairsCounted);
					gTExFeatureCounter.logInfo();
				}
				
			}
		}
		
		logger.info("Finished processing reads");
		gTExFeatureCounter.logInfo();
		
		double numberOfReads = gTExFeatureCounter.getTotalCount();
		List<FeatureCount> counts = gTExFeatureCounter.getCounts();
		List<Feature> newFeatures = new LinkedList<>();
		
		//Update features
		logger.info("Updating exons to include read counts and RPKM");
		for(FeatureCount fc : counts) {
			Feature f = fc.getFeature();
			
			Map<String,String> newAtts = new HashMap<>();
			newAtts.put("id", fc.getId());
			newAtts.put("reads", Double.toString(fc.getCount()));
			newAtts.put("RPKM", Double.toString(fc.getRPKM(numberOfReads)));
			
			Feature newF = GTFFeatureBuilder.addAttributesToFeature(f, newAtts);
			newFeatures.add(newF);
		}
		
		//Get Gene Counts
		
		logger.info("Writing GTF file: " + fileOut);
		try(GTFWriter writer = GTFWriter.getGTFWriterForFile(new File(fileOut))) {
			writer.write(newFeatures);
		}
	}
	
	
}
