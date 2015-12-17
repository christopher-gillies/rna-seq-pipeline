package org.kidneyomics.gtf;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.Map;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;

public class DefaultReadLogger implements ReadLogger {

	
	private File file;
	private BufferedWriter writer;
	private boolean isOpen = false;
	private Logger logger = LoggerFactory.getLogger(DefaultReadLogger.class);
	@Override
	public void setFile(File file) {
		if(file == null) {
			throw new IllegalArgumentException("file cannot be null");
		}
		this.file = file;
	}

	private void open() throws IOException {
		try {
			if(this.file == null) {
				throw new IllegalStateException("file cannot be null");
			}
			Charset charset = Charset.forName("US-ASCII");
			this.writer = Files.newBufferedWriter(file.toPath(), charset);
			isOpen = true;
		} catch(Exception e) {
			throw new IllegalStateException("Could not open writer for file");
		}
	}
	
	@Override
	public void logRead(String type, SAMRecord record, Map<Feature,Integer> mappedFeatures) {
		try {
			if(!isOpen) {
				open();
				//write header
				writer.write("#TYPE");
				writer.write("\t");
				writer.write("READ_NAME");
				writer.write("\t");
				writer.write("REFERENCE_NAME");
				writer.write("\t");
				writer.write("ALN_START");
				writer.write("\t");
				writer.write("CIGAR");
				writer.write("\t");
				writer.write("ALN_END");
				writer.write("\t");
				writer.write("READ_LENGTH");
				writer.write("\t");
				writer.write("MAPPED_FEATURES");
				writer.write("\t");
				writer.write("TAGS");
				writer.write("\n");
			} 
			
			writer.write(type);
			writer.write("\t");
			writer.write(record.getReadName());
			writer.write("\t");
			writer.write(record.getReferenceName());
			writer.write("\t");
			writer.write( "" + record.getAlignmentStart());
			writer.write("\t");
			writer.write(record.getCigarString());
			writer.write("\t");
			writer.write("" + record.getAlignmentEnd());
			writer.write("\t");
			writer.write("" + record.getReadLength());
			writer.write("\t");
			if(mappedFeatures == null || mappedFeatures.size() == 0) {
				writer.write(".");
			} else {
				for(Map.Entry<Feature, Integer> entry : mappedFeatures.entrySet()) {
					writer.write(entry.getKey().getAttribute("id"));
					writer.write("_BASES_");
					writer.write("" + entry.getValue() + ",");
				}
			}
			writer.write("\t");
			for(SAMTagAndValue tag : record.getAttributes()) {
				writer.write(tag.tag + ":" + tag.value.toString());
				writer.write("_");
			}
			writer.write("\n");
		} catch(Exception exception) {
			logger.error("Could not log read");
			exception.printStackTrace();
		}
	}

	@Override
	public void close() {
		try {
			if(isOpen) {
				writer.close();
				isOpen = false;
			}
		} catch(Exception e) {
			logger.error("Could not close logger");
		}
	}

}
