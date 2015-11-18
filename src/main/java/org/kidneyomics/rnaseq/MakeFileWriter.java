package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.stringtemplate.v4.ST;

@Component
public class MakeFileWriter {

	@Autowired
	ApplicationOptions applicationOptions;
	
	Logger logger;
	
	@Autowired
	MakeFileWriter(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	public void writeMakeFile(Mode mode) throws IOException {
		switch(mode) {
		case ALIGN:
			logger.info("Writing alignment commands");
			writeAlignCommands();
			break;
		case ERROR:
			break;
		}
	}
	
	private void writeAlignCommands() throws IOException {
		
		//Read samples
		
		//Create MakeFile
		
		MakeFile make = new MakeFile();
		
		String dirBase = applicationOptions.getOutputDirectory();
		
		File baseDir = new File(dirBase);
		
		if(!baseDir.exists()) {
			logger.info("Creating " + baseDir);
			baseDir.mkdirs();
		}
		
		//Generating Genome Index for 1st pass
		
		MakeEntry genomeIndex1stPassEntry = new MakeEntry();
		
		ST genomeIndex1stPass = new ST("<star> --runMode genomeGenerate --genomeDir <genomeDir> --genomeFastaFiles <reference> --runThreadN <n> --sjdbGTFfile <gtf>");
		genomeIndex1stPass.add("star", applicationOptions.getStar());
		
		//Make directory
		String outDir = dirBase + "/" + "genomeDirPassOne";
		File outDirRef = new File(outDir);
		if(!outDirRef.exists()) {
			logger.info("Creating " + outDir);
			outDirRef.mkdirs();
		}
		genomeIndex1stPass.add("genomeDir", outDir);
		
		genomeIndex1stPass.add("reference", applicationOptions.getReferenceSequence());
		
		genomeIndex1stPass.add("n", applicationOptions.getNumThreadsGenomeIndex());
		
		genomeIndex1stPass.add("gtf", applicationOptions.getGtf());
		
		genomeIndex1stPassEntry.setTarget("FIRST_PASS_GENOME_INDEX.OK");
		
		genomeIndex1stPassEntry.addCommand(genomeIndex1stPass.render());
		
		genomeIndex1stPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(genomeIndex1stPassEntry);
		
		
		
		
		
		//Write makefile
		String makefileText = make.toString();
		
		logger.info("Writing Makefile");
		
		FileUtils.write(new File(baseDir + "/Makefile"), makefileText);
	}
	
}
