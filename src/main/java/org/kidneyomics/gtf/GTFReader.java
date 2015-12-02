package org.kidneyomics.gtf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.util.StringUtils;

import htsjdk.samtools.util.CloseableIterator;

public class GTFReader implements Iterable<Feature>, Iterator<Feature>, CloseableIterator<Feature> {

	private BufferedReader reader;
	private String currentLine = null;
	private String nextLine = null;
	private Logger logger;
	
	private GTFReader(BufferedReader buffered) throws IOException { 
		
		this.reader = buffered;
		this.logger = LoggerFactory.getLogger(GTFReader.class);
		//Skip header lines
		while((this.nextLine = reader.readLine()) != null && this.nextLine.startsWith("#"));
	}
	

	public static GTFReader getGTFByFile(File f) throws IOException {
		
		
		BufferedReader buffered;
		if(f.getAbsolutePath().endsWith(".gz")) {
			InputStream fileStream = new FileInputStream(f);
			InputStream gzipStream = new GZIPInputStream(fileStream);
			Reader decoder = new InputStreamReader(gzipStream);
			buffered = new BufferedReader(decoder);
		} else {
			InputStream fileStream = new FileInputStream(f);
			Reader decoder = new InputStreamReader(fileStream);
			buffered = new BufferedReader(decoder);
		}
		
		GTFReader gtfReader = new GTFReader(buffered);
	
		return gtfReader;
	}
	
	@Override
	public Iterator<Feature> iterator() {
		return this;
	}

	@Override
	public boolean hasNext() {
		if(!StringUtils.isEmpty(nextLine)) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public Feature next() {
		this.currentLine = this.nextLine;
		try {
			//Returns: A String containing the contents of the line, not including any line-termination characters, or null if the end of the stream has been reached
			this.nextLine = reader.readLine();
		} catch (IOException e) {
			logger.error(e.getMessage());
		}
		if(StringUtils.isEmpty(nextLine)) {
			try {
				reader.close();
			} catch (IOException e) {
				logger.error(e.getMessage());
			}
		}
		//return null if the current line is null
		if(StringUtils.isEmpty(this.currentLine)) {
			return null;
		} else {
			return GTFFeatureBuilder.createFromLine(currentLine);
		}
	}

	@Override
	public void remove() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void close() {
		try {
			if(reader != null) {
				reader.close();
			}
		} catch (IOException e) {
			logger.error(e.getMessage());
		}
		
	}
	
	public List<Feature> readAllLines() {
		List<Feature> list = new LinkedList<Feature>();
		for(Feature f : this) {
			list.add(f);
		}
		return list;
		
	}

}
