package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ReadCountStatMerger {

	private ApplicationOptions applicationOptions;
	
	private Logger logger;
	
	@Autowired
	public ReadCountStatMerger(LoggerService loggerService, ApplicationOptions applicationOptions) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
	}
	
	
	public void mergeStatFiles() throws IOException {
		String statList = applicationOptions.getFileIn();
		String outfile = applicationOptions.getFileOut();
		
		List<String> lines = FileUtils.readLines(new File(statList));
		List<StatFile> statFiles = new LinkedList<>();
		for(String line : lines) {
			if(line.length() == 0) {
				continue;
			}
			String cols[] = line.split("\t");
			String id = cols[0];
			File file = new File(cols[1]);
			statFiles.add( new StatFile(file, id));
		}
		
		mergeStatFiles(new File(outfile), statFiles);
	}
	
	
	private void mergeStatFiles(File out, List<StatFile> statFiles) throws IOException {
		assert(statFiles.size() > 0);
		Collections.sort(statFiles);
		
		StatFile first = statFiles.get(0);
		
		StringBuilder sb = new StringBuilder();
		String[] cols = first.cols;
		String header = StringUtils.join(first.cols, "\t").replaceAll("#", "");
		
		sb.append("#ID");
		sb.append("\t");
		sb.append(header);
		sb.append("\n");
		for(StatFile sf : statFiles) {
			List<String> vals = new LinkedList<>();
			for(String col : cols) {
				vals.add(sf.data.get(col));
			}
			String line = StringUtils.join(vals,"\t");
			
			sb.append(sf.id);
			sb.append("\t");
			sb.append(line);
			sb.append("\n");
		}
		
		FileUtils.write(out, sb.toString());
		
	}
	
	private class StatFile implements Comparable<StatFile> {
		
		final String id;
		final String cols[];
		final String vals[];
		final HashMap<String,String> data = new HashMap<>();
		
		StatFile(File in, String id) throws IOException {
			this.id = id;
			
			List<String> lines = FileUtils.readLines(in);
			if(lines.size() < 2) {
				throw new RuntimeException(in.getAbsolutePath() + " is not formatted correctly");
			}
			
			cols = lines.get(0).split("\t");
			vals = lines.get(1).split("\t");
			
			if(cols.length != vals.length) {
				throw new RuntimeException(in.getAbsolutePath() + " is not formatted correctly. There are more columns than values");
			}
			
			for(int i = 0; i < cols.length; i++) {
				data.put(cols[i], vals[i]);
			}
		}

		@Override
		public int compareTo(StatFile o) {
			return this.id.compareTo(o.id);
		}
		
	}
}
