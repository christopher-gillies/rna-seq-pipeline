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

public class GTFReader implements Iterable<Feature>, Iterator<Feature>, CloseableIterator<Feature>, AutoCloseable {

	private BufferedReader reader;
	private Feature currentFeature = null;
	private String nextLine = null;
	private Feature nextFeature = null;
	private Logger logger;
	private boolean noVersion = false;
	private List<GTFFeatureFilter> filters = new LinkedList<GTFFeatureFilter>();
	private long numFilteredOut = 0;
	
	private GTFReader(BufferedReader buffered, boolean noEnsembleVersion) throws IOException { 
		this.noVersion = noEnsembleVersion;
		this.reader = buffered;
		this.logger = LoggerFactory.getLogger(GTFReader.class);
		//Skip header lines
		while((this.nextLine = reader.readLine()) != null && this.nextLine.startsWith("#"));
		//no filters added at this point
		getNextFeature();
	}
	
	private Feature getNextFeature() {
		nextFeature = null;
		//next line is pointing to the next element
		if(!StringUtils.isEmpty(nextLine)) {
			if(this.noVersion) {
				nextFeature = GTFFeatureBuilder.createFromLine(nextLine,true);
			} else {
				nextFeature = GTFFeatureBuilder.createFromLine(nextLine);
			}
			
			//read next line
			try {
				//Returns: A String containing the contents of the line, not including any line-termination characters, or null if the end of the stream has been reached
				this.nextLine = reader.readLine();
			} catch (IOException e) {
				logger.error(e.getMessage());
			}
			
			boolean keep = keep(nextFeature);
			
			//if we are not keeping this feature then get the next feature 
			if(!keep) {
				numFilteredOut++;
				getNextFeature();
			}
		}
		return nextFeature;
		
	}

	private boolean keep(Feature f) {
		boolean keep = true;
		for(GTFFeatureFilter filter : filters) {
			boolean filterRes = filter.filter(f);
			keep = keep && filterRes;
		}
		return keep;
	}
	private static BufferedReader getBufferedReader(File f) throws IOException {
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
		return buffered;
	}
	
	public static GTFReader getGTFByFile(File f) throws IOException {
		
		
		BufferedReader buffered = getBufferedReader(f);
		GTFReader gtfReader = new GTFReader(buffered,false);
	
		return gtfReader;
	}
	
	
	public static GTFReader getGTFByFileNoEnsemblVersion(File f) throws IOException {
		
		
		BufferedReader buffered = getBufferedReader(f);
		
		GTFReader gtfReader = new GTFReader(buffered,true);
	
		return gtfReader;
	}
	
	@Override
	public Iterator<Feature> iterator() {
		return this;
	}

	@Override
	public boolean hasNext() {
		if(this.nextFeature != null) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public Feature next() {
		this.currentFeature = this.nextFeature;
		getNextFeature();
		return currentFeature;
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
	
	public GTFReader addFilter(GTFFeatureFilter filter) {
		this.filters.add(filter);
		/*
		 * After adding a filter get the next feature meeting it
		 */
		if(!keep(this.nextFeature)) {
			numFilteredOut++;
			getNextFeature();
		}
		return this;
	}

	public long getNumberFilteredOut() {
		return this.numFilteredOut;
	}
}
