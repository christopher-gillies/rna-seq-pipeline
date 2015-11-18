package org.kidneyomics.rnaseq;

import org.apache.commons.cli.ParseException;

public interface OptionProcessor {
	void processInputs(String[] args) throws ParseException;
}
