package org.kidneyomics.rnaseq;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;

public class QuantificationUtils {

	public static void writeQuantificationMatrix(Collection<? extends Quantification> listOfQuantifications, String outmatrix) throws IOException {
		Path p = Paths.get(outmatrix);
		BufferedWriter bf = Files.newBufferedWriter(p,Charset.defaultCharset());
		
		int index = 0;
		for(Quantification tq : listOfQuantifications) {
			if(index == 0) {
				bf.append(tq.printHeader());
				bf.append("\n");
			}
			index++;
			tq.appendTo(bf);
			bf.append("\n");
		}
		
		bf.close();
	}
}
