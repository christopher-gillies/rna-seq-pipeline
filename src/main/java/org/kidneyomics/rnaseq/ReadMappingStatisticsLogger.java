package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;

public class ReadMappingStatisticsLogger {
	
	public static void writeStats(File file, FeatureCounter fc) throws IOException {
		double ambiguousReadCount = fc.getAmbiguousReadCount();
		double mappedReadCount = fc.getMappedReadCount();
		double totalReads = fc.getTotalCount();
		long partiallyMappedReads = fc.getNumberOfPartiallyUnmappedReads();
		double unmappedReadCount = fc.getUnmappedReadCount();
		
		StringBuilder sb = new StringBuilder();
		
		
		sb.append("#TOTAL_READS");
		sb.append("\t");
		sb.append("AMBIGUOUS_READS");
		sb.append("\t");
		sb.append("AMBIGUOUS_FRACTION");
		sb.append("\t");
		sb.append("UNMAPPED_READS");
		sb.append("\t");
		sb.append("UNMAPPED_FRACTION");
		sb.append("\t");
		sb.append("PARTIALLY_MAPPED_READS");
		sb.append("\t");
		sb.append("PARTIALLY_MAPPED_FRACTION");
		sb.append("\t");
		sb.append("MAPPED_READS");
		sb.append("\t");
		sb.append("MAPPED_FRACTION");
		sb.append("\n");
		
		
		
		sb.append(totalReads);
		sb.append("\t");
		sb.append(ambiguousReadCount);
		sb.append("\t");
		sb.append(ambiguousReadCount / totalReads);
		sb.append("\t");
		sb.append(unmappedReadCount);
		sb.append("\t");
		sb.append(unmappedReadCount / totalReads);
		sb.append("\t");
		sb.append(partiallyMappedReads);
		sb.append("\t");
		sb.append(partiallyMappedReads / totalReads);
		sb.append("\t");
		sb.append(mappedReadCount);
		sb.append("\t");
		sb.append(mappedReadCount / totalReads);
		sb.append("\n");
		
		FileUtils.write(file, sb.toString());
	}
	
	
}
