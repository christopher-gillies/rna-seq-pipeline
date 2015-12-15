package org.kidneyomics.gtf;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;

import org.biojava.nbio.genome.parsers.gff.Feature;

public class GTFWriter implements AutoCloseable {
	
	BufferedWriter bf;
	
	private GTFWriter(File f) throws IOException {
	
		Path p = f.toPath();
		bf = Files.newBufferedWriter(p,Charset.defaultCharset());
		
	}
	
	public void write(Collection<Feature> features) throws IOException {
		for(Feature feature : features) {
			write(feature);
		}
	}
	
	public void write(Feature feature) throws IOException {
		bf.write(GTFFeatureRenderer.render(feature));
		bf.write("\n");
	}
	
	
	
	
	public static GTFWriter getGTFWriterForFile(File f) throws IOException {
		return new GTFWriter(f);
	}

	@Override
	public void close() throws IOException {
		bf.close();
	}
}
