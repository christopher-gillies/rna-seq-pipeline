package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import org.apache.commons.io.FileUtils;

public class ReadMappingStatisticsLogger {
	
	private static final DecimalFormat df = new DecimalFormat("#.#");
	static  {
		df.setMaximumFractionDigits(100);
	}
	
	public static void writeStats(File file, FeatureCounter fc) throws IOException {
		double ambiguousReadCount = fc.getAmbiguousReadCount();
		double mappedReadCount = fc.getMappedReadCount();
		double totalReads = fc.getTotalCount();
		long partiallyMappedReads = fc.getNumberOfPartiallyUnmappedReads();
		long filtered = fc.getNumberOfFilteredReads();
		double unmappedReadCount = fc.getUnmappedReadCount();
		
		StringBuilder sb = new StringBuilder();
		
		
		sb.append("#TOTAL_READS");
		sb.append("\t");
		sb.append("AMBIGUOUS_READS");
		sb.append("\t");
		sb.append("AMBIGUOUS_FRACTION");
		sb.append("\t");
		sb.append("FILTERED_OUT_READS");
		sb.append("\t");
		sb.append("FILTERED_OUT_FRACTION");
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
		
		
		
		sb.append( df.format(totalReads) );
		sb.append("\t");
		sb.append( df.format(ambiguousReadCount));
		sb.append("\t");
		sb.append( df.format(ambiguousReadCount / totalReads) );
		sb.append("\t");
		sb.append( df.format(filtered));
		sb.append("\t");
		sb.append( df.format(filtered / totalReads) );
		sb.append("\t");
		sb.append( df.format( unmappedReadCount ));
		sb.append("\t");
		sb.append( df.format( unmappedReadCount / totalReads));
		sb.append("\t");
		sb.append(partiallyMappedReads);
		sb.append("\t");
		sb.append(df.format( partiallyMappedReads / totalReads));
		sb.append("\t");
		sb.append( df.format( mappedReadCount));
		sb.append("\t");
		sb.append(df.format(mappedReadCount / totalReads));
		sb.append("\n");
		
		FileUtils.write(file, sb.toString());
	}
	
	
}
