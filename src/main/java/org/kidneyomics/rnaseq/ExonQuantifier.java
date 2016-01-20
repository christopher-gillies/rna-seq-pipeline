package org.kidneyomics.rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.kidneyomics.gtf.DefaultReadLogger;
import org.kidneyomics.gtf.FeatureCount;
import org.kidneyomics.gtf.GTFFeatureBuilder;
import org.kidneyomics.gtf.GTFFeatureUtil;
import org.kidneyomics.gtf.GTFWriter;
import org.kidneyomics.gtf.ReadLogger;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

@Component
public class ExonQuantifier implements ApplicationCommand {

	
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
		
		String readLogFile = applicationOptions.getReadLogFile();
		if(readLogFile != null) {
			ReadLogger readLogger = new DefaultReadLogger();
			readLogger.setFile(new File(readLogFile));
			gTExFeatureCounter.setReadLogger(readLogger);
		}
		
		//Read GTF file
		//FindOverlappingFeatures findOverlappingFeatures = new FindOverlappingFeatures();
		//FeatureCounter featureCounter = new GTExFeatureCounter(findOverlappingFeatures, loggerService)
		
		gTExFeatureCounter.addFilter(new MaxAlignmentDistanceFilter(applicationOptions.getMaxEditDistance()));
		
		gTExFeatureCounter.buildFeatures(new File(gtfFile), "exon");
		
		
		//Process BAM
		File fin = new File(fileIn);
		int pairsCounted = 0;
		try(BAMProcessor processor = BAMProcessorImplDict.getBAMProcessor(fin)) {
			
			SAMRecordPair pair = null;
			while( ( pair = processor.getNextReadPair()) != null) {
				gTExFeatureCounter.count(pair);
				pairsCounted++;
				if(pairsCounted % 100000 == 0) {
					logger.info("Read pairs counted: " + pairsCounted);
					gTExFeatureCounter.logInfo();
				}
				
			}
		}
		
		logger.info("Finished processing reads");
		gTExFeatureCounter.logInfo();
		
		gTExFeatureCounter.getReadLogger().close();
		
		logger.info("Updating exons to include read counts and RPKM");
		List<Feature> exonCounts = gTExFeatureCounter.getCounts();
		List<Feature> geneCounts =  gTExFeatureCounter.getGeneCounts();
		List<Feature> allCounts = new ArrayList<>(exonCounts.size() + geneCounts.size());
		allCounts.addAll(exonCounts);
		allCounts.addAll(geneCounts);
		//Get Gene Counts
		
		//Write new features
		logger.info("Sorting exon and gene results by chromosome position");
		GTFFeatureUtil.sortFeatures(allCounts);
	
		logger.info("Writing GTF file: " + fileOut);
		try(GTFWriter writer = GTFWriter.getGTFWriterForFile(new File(fileOut))) {
			writer.write(allCounts);
		}
		
		logger.info("Writing log files");
		//Log statistics
		ReadMappingStatisticsLogger.writeStats(new File(fileOut + ".stats"), gTExFeatureCounter);
		
		logger.info("Write out removed features");
		writeOutRemovedExons(gTExFeatureCounter.getRemovedFeatures(),fileOut);
	}
	
	//write out exons that overlap
	private void writeOutRemovedExons(Collection<Feature> removedFeatures, String fileOut) throws IOException {
		
		try(BufferedWriter writer = Files.newBufferedWriter(Paths.get(fileOut + ".exons.removed"), Charset.defaultCharset(), StandardOpenOption.CREATE)) {
		
			for(Feature feature : removedFeatures) {
				String gene = feature.getAttribute("gene_id");
				String transcript = feature.getAttribute("transcript_id");
				String chr = feature.seqname();
				Integer start = feature.location().bioStart();
				Integer end = feature.location().bioEnd();
				
				writer.append(transcript);
				writer.append("\t");
				
				writer.append(gene);
				writer.append("\t");
				
				writer.append(chr);
				writer.append("\t");
				
				writer.append(start.toString());
				writer.append("\t");
				
				writer.append(end.toString());
				
				writer.append("\n");
			}
		}
	}


	@Override
	public void doWork() throws Exception {
		quantify();
	}
	
	
}
