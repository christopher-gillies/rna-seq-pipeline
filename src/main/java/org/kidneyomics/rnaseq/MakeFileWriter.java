package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;

import org.apache.commons.io.FileUtils;
import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;
import org.springframework.util.StringUtils;
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
	
	public void writeMakeFile(Mode mode) throws Exception {
		switch(mode) {
		case ALIGN:
			logger.info("Writing alignment commands");
			writeAlignCommands();
			break;
		case ERROR:
			break;
		}
	}
	
	private void writeAlignCommands() throws Exception {
		
		//Read samples
		
		//Create MakeFile
		
		MakeFile make = new MakeFile();
		
		String dirBase = applicationOptions.getOutputDirectory();
		
		File baseDir = new File(dirBase);
		
		if(!baseDir.exists()) {
			logger.info("Creating " + baseDir);
			baseDir.mkdirs();
		}
		
		/*
		 * 
		 * 
		 * Generating Genome Index for 1st pass
		 * 
		 * 
		 * 
		 */
		
		MakeEntry genomeIndex1stPassEntry = new MakeEntry();
		
		
		genomeIndex1stPassEntry.setComment("Generate the genome for first pass alignment");
		
		ST genomeIndex1stPass = new ST("<star> --runMode genomeGenerate --genomeDir <genomeDir> --genomeFastaFiles <reference> --runThreadN <n> --sjdbGTFfile <gtf> --sjdbOverhang <overhang>");
		genomeIndex1stPass.add("star", applicationOptions.getStar());
		
		//Make directory
		String genomeOutDir = dirBase + "/" + "genomeDirPassOne";
		File outDirRef = new File(genomeOutDir);
		if(!outDirRef.exists()) {
			logger.info("Creating " + genomeOutDir);
			outDirRef.mkdirs();
		}
		genomeIndex1stPass.add("genomeDir", genomeOutDir);
		
		genomeIndex1stPass.add("reference", applicationOptions.getReferenceSequence());
		
		genomeIndex1stPass.add("n", applicationOptions.getNumThreadsGenomeIndex());
		
		genomeIndex1stPass.add("gtf", applicationOptions.getGtf());
		
		genomeIndex1stPass.add("overhang", applicationOptions.getSjdbOverhang());
		
		genomeIndex1stPassEntry.setTarget(dirBase + "/FIRST_PASS_GENOME_INDEX.OK");
		
		genomeIndex1stPassEntry.addCommand(genomeIndex1stPass.render());
		
		genomeIndex1stPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(genomeIndex1stPassEntry);
		
		
		/*
		 * 
		 * 
		 * Load genome command
		 * 
		 * 
		 * 
		 */
		MakeEntry loadGenomeEntry = new MakeEntry();
		loadGenomeEntry.setComment("Load the genome into RAM");
		loadGenomeEntry.addDependency(genomeIndex1stPassEntry);
		loadGenomeEntry.setTarget(baseDir + "/GENOME_LOADED.OK");
		
		ST loadGenome = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndExit");
		loadGenome.add("star", applicationOptions.getStar());
		loadGenome.add("genomeDir", genomeOutDir);
		
		
		loadGenomeEntry.addCommand(loadGenome.render());
		loadGenomeEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(loadGenomeEntry);
		
		
		/*
		 * 
		 * 
		 * Generating Alignment for first pass
		 * 
		 * 
		 * 
		 */
		
		//Read files
		Collection<Sample> samples = Sample.getFastqFileList(new File(applicationOptions.getFastqFiles()));
		
		
		
		Collection<MakeEntry> pass1Dependencies = new LinkedList<MakeEntry>();
		
	
		Collection<String> sjdbFiles = new LinkedList<String>();
		
		for(Sample sample : samples) {
			MakeEntry firstPassAlignEntry = new MakeEntry();
			firstPassAlignEntry.setComment("Load the genome into RAM");
			firstPassAlignEntry.addDependency(genomeIndex1stPassEntry);
			firstPassAlignEntry.addDependency(loadGenomeEntry);
			firstPassAlignEntry.setTarget(dirBase +  "/" + sample.getSampleId()  + ".FIRST_PASS.OK");
			
			String sampleDir = dirBase + "/" + sample.getSampleId() + "_1/";
			
			ST firstPassAlign = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndKeep --readFilesIn <files> --readFilesCommand <uncompress> --outFileNamePrefix <outdir>");
			firstPassAlign.add("star", applicationOptions.getStar());
			firstPassAlign.add("genomeDir", genomeOutDir);
			firstPassAlign.add("files", StringUtils.collectionToDelimitedString(sample.getFastqFiles(), " "));
			firstPassAlign.add("uncompress", applicationOptions.getUncompressCommand());
			firstPassAlign.add("outdir", sampleDir);
			
			
			String sjdbFile = sampleDir + "/SJ.out.tab";
			sjdbFiles.add(sjdbFile);
			
			firstPassAlignEntry.addCommand("mkdir -p " + sampleDir);
			firstPassAlignEntry.addCommand(firstPassAlign.render());
			firstPassAlignEntry.addCommand("touch $@");
			
			// add to make file
			make.addMakeEntry(firstPassAlignEntry);
			
			//add to dependencies
			pass1Dependencies.add(firstPassAlignEntry);
		}
		
		
		/*
		 * 
		 * 
		 * Generating Alignment for second pass
		 * 
		 * 
		 * 
		 */
		
		Collection<MakeEntry> pass2Dependencies = new LinkedList<MakeEntry>();
		
		for(Sample sample : samples) {
			MakeEntry secondPassAlignEntry = new MakeEntry();
			secondPassAlignEntry.setComment("Load the genome into RAM");
			secondPassAlignEntry.addDependency(genomeIndex1stPassEntry);
			secondPassAlignEntry.addDependency(loadGenomeEntry);
			secondPassAlignEntry.addDependencies(pass1Dependencies);
			
			secondPassAlignEntry.setTarget(dirBase +  "/" + sample.getSampleId()  + ".SECOND_PASS.OK");
			
			String sampleDir = dirBase + "/" + sample.getSampleId() + "/";
			
			ST SecondPassAlign = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndKeep --readFilesIn <files> --readFilesCommand <uncompress> --sjdbFileChrStartEnd <sjdbs> --outFileNamePrefix <outdir>");
			SecondPassAlign.add("star", applicationOptions.getStar());
			SecondPassAlign.add("genomeDir", genomeOutDir);
			SecondPassAlign.add("files", StringUtils.collectionToDelimitedString(sample.getFastqFiles(), " "));
			SecondPassAlign.add("sjdbs", StringUtils.collectionToDelimitedString(sjdbFiles, " "));
			SecondPassAlign.add("uncompress", applicationOptions.getUncompressCommand());
			SecondPassAlign.add("outdir", sampleDir);
			
			
			
			
			secondPassAlignEntry.addCommand("mkdir -p" + sampleDir);
			secondPassAlignEntry.addCommand(SecondPassAlign.render());
			secondPassAlignEntry.addCommand("touch $@");
			
			// add to make file
			make.addMakeEntry(secondPassAlignEntry);
			
			//add to dependencies
			pass2Dependencies.add(secondPassAlignEntry);
		}
		
		
		//Should be done last
		/*
		 * 
		 * 
		 * Remove genome command
		 * 
		 * 
		 * 
		 */
		MakeEntry removeGenomeEntry = new MakeEntry();
		removeGenomeEntry.setComment("Remove the genome from RAM");
		removeGenomeEntry.addDependency(genomeIndex1stPassEntry);
		removeGenomeEntry.addDependency(loadGenomeEntry);
		
		removeGenomeEntry.addDependencies(pass1Dependencies);
		removeGenomeEntry.addDependencies(pass2Dependencies);
		
		removeGenomeEntry.setTarget(baseDir + "/GENOME_REMOVED.OK");
		
		ST removeGenome = new ST("<star> --genomeDir <genomeDir> --genomeLoad Remove");
		removeGenome.add("star", applicationOptions.getStar());
		removeGenome.add("genomeDir", genomeOutDir);
		
		
		removeGenomeEntry.addCommand(removeGenome.render());
		removeGenomeEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(removeGenomeEntry);
		
		
		
		/*
		 * 
		 * 
		 * Write makefile
		 * 
		 * 
		 */
		String makefileText = make.toString();
		
		logger.info("Writing Makefile");
		
		FileUtils.write(new File(baseDir + "/Makefile"), makefileText);
	}
	
}
