package org.kidneyomics.rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class MapIdentifiersForExpressionMatrix implements ApplicationCommand {

	
	Logger logger;
	ApplicationOptions applicationOptions;
	
	@Autowired
	public MapIdentifiersForExpressionMatrix(LoggerService loggerService, ApplicationOptions applicationOptions) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
	}
	
	@Override
	public void doWork() throws Exception {
		
		//lsat column appears to be strand for gene, transcript, and exon
		final String expressionMatrix = applicationOptions.getExpressionMatrix();
		
		//OLD_ID\tNEW_ID
		final String infile = applicationOptions.getFileIn();
		final String outfile = applicationOptions.getFileOut();
		final String lastCol = applicationOptions.getLastNonSampleColInExpressionMatrix();
		
		//Map
		HashMap<String,String> map = new HashMap<>();
		
		//Read infile
		List<String> ids = FileUtils.readLines(new File(infile));
		for(String line : ids) {
			if( StringUtils.isEmpty(line)) {
				continue;
			}
			String[] cols = line.split("\t");
			
			if(cols.length != 2) {
				throw new IllegalArgumentException(infile + " should have two colums OLD_ID\tNEW_ID");
			}
			
			map.put(cols[0], cols[1]);
			
		}
		
		//we need to read the expression matrix and specify a new header
		//test file is ./src/test/resources/test_genes.txt  
		
		//autocloseable
		try(BufferedReader reader = Files.newBufferedReader(Paths.get(expressionMatrix), Charset.defaultCharset())) {
			
			try(BufferedWriter writer = Files.newBufferedWriter(Paths.get(outfile), Charset.defaultCharset(), StandardOpenOption.CREATE)) {
			
				String header = null;
				header = reader.readLine();
				String[] headerCols = header.split("\t");
				
				List<String> newHeader = new ArrayList<>(headerCols.length);
				
				//find the index of the last column
				final int lastColIndex = ArrayUtils.indexOf(headerCols, lastCol);
				
				//save all the items up to and including the alst index in the new header
				for(int i = 0; i <= lastColIndex; i++) {
					newHeader.add( headerCols[i]  );
				}
				
				//map the new items into the 
				for(int i = lastColIndex + 1; i < headerCols.length; i++) {
					//map columns
					String oldId = headerCols[i];
					if(map.containsKey(oldId)) {
						newHeader.add( map.get(oldId) );
					} else {
						newHeader.add(oldId);
					}
				}
				
				//write the new header
				
				String newHeaderLine = StringUtils.join(newHeader,"\t");
				writer.write(newHeaderLine);
				writer.write("\n");
				
				//write the rest of the file
				//the end of input is -1
				int nextChar = -1;
				while( (nextChar = reader.read()) != -1) {
					writer.write(nextChar);
				}
			}
		
		}
		
	}

}
