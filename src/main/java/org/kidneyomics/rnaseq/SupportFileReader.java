package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.List;

public interface SupportFileReader<T> {
	public List<T> readFile(File f) throws IOException;
	public List<String> getIds(List<T> items);
}
