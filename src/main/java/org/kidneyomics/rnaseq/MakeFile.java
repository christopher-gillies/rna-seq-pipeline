package org.kidneyomics.rnaseq;

import java.util.Collection;
import java.util.LinkedList;

import org.apache.log4j.Logger;
import org.stringtemplate.v4.ST;
import org.stringtemplate.v4.STGroup;
import org.stringtemplate.v4.STGroupFile;

public class MakeFile implements TemplateStringWriter {
	
	private Collection <MakeEntry> makeEntries;
	private final String makefileTemplate = "st/makefile.stg";
	
	Logger logger = Logger.getLogger(MakeFile.class);
	
	public MakeFile() {
		this.makeEntries = new LinkedList<MakeEntry>();		
		logger.info("MakeFile created");
		logger.info(makefileTemplate);
	}
	
	public Collection<MakeEntry> getMakeEntries() {
		return makeEntries;
	}

	public MakeFile addMakeEntry(MakeEntry entry) {
		this.makeEntries.add(entry);
		return this;
	}

	public String writeTemplateToString() {
		
		STGroup group = new STGroupFile(makefileTemplate);
		ST st = group.getInstanceOf("makefile");
		
		st.add("entries", this.makeEntries);
		
		return st.render();
	}
	
	@Override
	public String toString() {
		return writeTemplateToString();
	}

	
	
	
}
