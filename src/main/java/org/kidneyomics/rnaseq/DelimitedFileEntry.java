package org.kidneyomics.rnaseq;

import java.io.BufferedWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;


class DelimitedFileEntry {
	
	private DelimitedFileEntry(String sampleId, String[] cols) {
		map = new HashMap<>();
		this.sampleId = sampleId;
		this.cols = cols;
	};
	
	private final String sampleId;
	
	private final HashMap<String,String> map;
	
	private final String[] cols;
	
	static DelimitedFileEntry getDelimitedFileEntryFromDelimitedLine(String sampleId, String header, String values, String delimiter) {
		return DelimitedFileEntry.getDelimitedFileEntry(sampleId, header.split(delimiter), values.split(delimiter));
	}
	
	static DelimitedFileEntry getDelimitedFileEntry(String sampleId, String[] header, String[] values) {
		
		if(header.length != values.length) {
			throw new RuntimeException("There are a different number of header columns and values." + "\n" + StringUtils.join(header,",") + "\n" + StringUtils.join(values,",")) ;
		}
		
		
		DelimitedFileEntry d = new DelimitedFileEntry(sampleId, header);
		
		//store values in hash
		for(int i = 0; i < header.length; i++) {
			d.map.put(header[i], values[i]);
		}
		
		return d;
	}
	
	/**
	 * 
	 * @param entries
	 * @return union of columns across all entries. trying best to preserve order
	 */
	static List<String> getAllKeys(List<DelimitedFileEntry> entries) {
		List<String> result = new LinkedList<String>();
		Set<String> tmp = new HashSet<>();
		//iterate through all entries and cols for each entry
		//add new column if we have not already added it
		for(DelimitedFileEntry entry : entries) {
			for(String col : entry.cols) {
				if(!tmp.contains(col)) {
					tmp.add(col);
					result.add(col);
				}
			}
		}
		
		
		return result;
	}
	
	/**
	 * 
	 * @param cols
	 * @param delimiter
	 * @return the cols separated by the delimiter
	 */
	static String headerToString(List<String> cols, String delimiter) {
		LinkedList<String> toJoin = new LinkedList<>();
		toJoin.add("SAMPLE_ID");
		toJoin.addAll(cols);
		return StringUtils.join(toJoin, delimiter);
	}
	
	/**
	 * 
	 * @param entries
	 * @param file
	 * @param delimited
	 * 
	 * writer the entries to the file with delimiter delimited
	 */
	static void writeToFile(List<DelimitedFileEntry> entries, String file, String delimited) {
		try(BufferedWriter writer = Files.newBufferedWriter(Paths.get(file), Charset.defaultCharset(), StandardOpenOption.CREATE)) {
			List<String> headerVals = DelimitedFileEntry.getAllKeys(entries);
			writer.write(  DelimitedFileEntry.headerToString(headerVals, delimited) );
			writer.write("\n");
			for(DelimitedFileEntry dfe : entries) {
				writer.write( dfe.toString(headerVals, delimited));
				writer.write("\n");
			}
		} catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @param cols
	 * @param delimiter
	 * @return the values for each col separated by the delimiter
	 */
	String toString(List<String> cols, String delimiter) {
		StringBuilder sb = new StringBuilder();
		Iterator<String> iter = cols.iterator();
		
		sb.append(this.getSampleId());
		if(iter.hasNext()) {
			sb.append(delimiter);
		}
		
		while(iter.hasNext()) {
			String col = iter.next();
			
			String val = this.getValue(col);
			
			sb.append(val);
			
			if(iter.hasNext()) {
				sb.append(delimiter);
			}
		}
		
		return sb.toString();
		
	}
	
	@Override
	public String toString() {
		
		return toString(Arrays.asList(this.cols),"\t");
	}
	
	String getSampleId() {
		return sampleId;
	}
	
	String getValue(String col) {
		String res = null;
		if(map.containsKey(col)) {
			res = map.get(col);
		} else {
			res = "NA";
		}
		return res;
	}
}
