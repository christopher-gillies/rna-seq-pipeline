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
		
		MakeEntry genomeIndexFirstPassEntry = new MakeEntry();

		genomeIndexFirstPassEntry.setComment("Generate the genome for first pass alignment");
		
		ST genomeIndexFirstPass = new ST("<star> --runMode genomeGenerate --genomeDir <genomeDir> --genomeFastaFiles <reference> --runThreadN <n> --sjdbGTFfile <gtf> --sjdbOverhang <overhang>");
		genomeIndexFirstPass.add("star", applicationOptions.getStar());
		
		//Make directory
		String genomeOutDir = dirBase + "/" + "genomeDirPassOne";
		File outDirRef = new File(genomeOutDir);
		if(!outDirRef.exists()) {
			logger.info("Creating " + genomeOutDir);
			outDirRef.mkdirs();
		}
		genomeIndexFirstPass.add("genomeDir", genomeOutDir);
		
		genomeIndexFirstPass.add("reference", applicationOptions.getReferenceSequence());
		
		genomeIndexFirstPass.add("n", applicationOptions.getNumThreadsGenomeIndex());
		
		genomeIndexFirstPass.add("gtf", applicationOptions.getGtf());
		
		genomeIndexFirstPass.add("overhang", applicationOptions.getSjdbOverhang());
		
		genomeIndexFirstPassEntry.setTarget(dirBase + "/FIRST_PASS_GENOME_INDEX.OK");
		
		genomeIndexFirstPassEntry.addCommand(genomeIndexFirstPass.render());
		
		genomeIndexFirstPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(genomeIndexFirstPassEntry);
		
		
		/*
		 * 
		 * 
		 * Load genome for first pass command
		 * 
		 * 
		 * 
		 */
		MakeEntry loadGenomeFirstPassEntry = new MakeEntry();
		
		loadGenomeFirstPassEntry.setComment("Load the genome into RAM");
		loadGenomeFirstPassEntry.addDependency(genomeIndexFirstPassEntry);
		loadGenomeFirstPassEntry.setTarget(baseDir + "/GENOME_LOADED_FIRST_PASS.OK");
		
		ST loadGenome = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndExit");
		loadGenome.add("star", applicationOptions.getStar());
		loadGenome.add("genomeDir", genomeOutDir);
		
		
		loadGenomeFirstPassEntry.addCommand(loadGenome.render());
		loadGenomeFirstPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(loadGenomeFirstPassEntry);
		
		
		/*
		 * 
		 * 
		 * Generating alignment for first pass
		 * 
		 * 
		 * 
		 */
		
		//Read files
		Collection<Sample> samples = Sample.getFastqFileList(new File(applicationOptions.getFastqFiles()));
		
		//Dependencies for pass one
		Collection<MakeEntry> pass1Dependencies = new LinkedList<MakeEntry>();
	
		//splice junction databases
		Collection<String> sjdbFiles = new LinkedList<String>();
		
		for(Sample sample : samples) {
			MakeEntry firstPassAlignEntry = new MakeEntry();
			firstPassAlignEntry.setComment("Run first pass alignment for sample " + sample.getSampleId());
			firstPassAlignEntry.addDependency(genomeIndexFirstPassEntry);
			firstPassAlignEntry.addDependency(loadGenomeFirstPassEntry);
			firstPassAlignEntry.setTarget(dirBase +  "/" + sample.getSampleId()  + ".FIRST_PASS.OK");
			
			String sampleDir = dirBase + "/" + sample.getSampleId() + "_1/";
			
			ST firstPassAlign = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndKeep --readFilesIn <files> --readFilesCommand <uncompress> --outFileNamePrefix <outdir> --outSJfilterCountUniqueMin 4 2 2 2 --outSJfilterCountTotalMin 4 2 2 2");
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
		 * Remove genome after first pass command 
		 * 
		 * 
		 * 
		 */
		MakeEntry removeGenomeFirstPassEntry = new MakeEntry();
		removeGenomeFirstPassEntry.setComment("Remove the genome from RAM");
		removeGenomeFirstPassEntry.addDependency(genomeIndexFirstPassEntry);
		removeGenomeFirstPassEntry.addDependency(loadGenomeFirstPassEntry);
		
		removeGenomeFirstPassEntry.addDependencies(pass1Dependencies);
		
		removeGenomeFirstPassEntry.setTarget(baseDir + "/GENOME_REMOVED_FIRST_PASS.OK");
		
		ST removeGenomeFirstPass = new ST("<star> --genomeDir <genomeDir> --genomeLoad Remove");
		removeGenomeFirstPass.add("star", applicationOptions.getStar());
		removeGenomeFirstPass.add("genomeDir", genomeOutDir);
		
		
		removeGenomeFirstPassEntry.addCommand(removeGenomeFirstPass.render());
		removeGenomeFirstPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(removeGenomeFirstPassEntry);
		
		
		
		/*
		 * 
		 * 
		 * Generating Genome Index for 2nd pass
		 * 
		 * 
		 * 
		 */
		
		MakeEntry genomeIndexSecondPassEntry = new MakeEntry();

		genomeIndexSecondPassEntry.setComment("Generate the genome for second pass alignment");
		
		
		genomeIndexSecondPassEntry.addDependencies(pass1Dependencies);
		genomeIndexSecondPassEntry.addDependency(removeGenomeFirstPassEntry);
		
		ST genomeIndexSecondPass = new ST("<star> --runMode genomeGenerate --genomeDir <genomeDir> --genomeFastaFiles <reference> --runThreadN <n> --sjdbGTFfile <gtf> --sjdbFileChrStartEnd <sjdbs> --sjdbOverhang <overhang>");
		genomeIndexSecondPass.add("star", applicationOptions.getStar());
		

		//Make directory
		String genomeOutDirPass2 = dirBase + "/" + "genomeDirPassTwo";
		File outDirRefPass2 = new File(genomeOutDirPass2);
		if(!outDirRefPass2.exists()) {
			logger.info("Creating " + genomeOutDirPass2);
			outDirRefPass2.mkdirs();
		}
		genomeIndexSecondPass.add("genomeDir", genomeOutDirPass2);
		
		genomeIndexSecondPass.add("reference", applicationOptions.getReferenceSequence());
		
		genomeIndexSecondPass.add("n", applicationOptions.getNumThreadsGenomeIndex());
		
		genomeIndexSecondPass.add("gtf", applicationOptions.getGtf());
		
		genomeIndexSecondPass.add("overhang", applicationOptions.getSjdbOverhang());
		
		genomeIndexSecondPass.add("sjdbs", StringUtils.collectionToDelimitedString(sjdbFiles, " "));
		
		genomeIndexSecondPassEntry.setTarget(dirBase + "/SECOND_PASS_GENOME_INDEX.OK");
		
		genomeIndexSecondPassEntry.addCommand(genomeIndexSecondPass.render());
		
		genomeIndexSecondPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(genomeIndexSecondPassEntry);
		
		
		/*
		 * 
		 * 
		 * Load genome for second pass command
		 * 
		 * 
		 * 
		 */
		MakeEntry loadGenomeSecondPassEntry = new MakeEntry();
		
		loadGenomeSecondPassEntry.setComment("Load the second pass genome into RAM");
		loadGenomeSecondPassEntry.addDependency(genomeIndexSecondPassEntry);
		loadGenomeSecondPassEntry.setTarget(baseDir + "/GENOME_LOADED_SECOND_PASS.OK");
		
		ST loadGenomeSecondPass = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndExit");
		loadGenomeSecondPass.add("star", applicationOptions.getStar());
		loadGenomeSecondPass.add("genomeDir", genomeOutDirPass2);
		
		
		loadGenomeSecondPassEntry.addCommand(loadGenomeSecondPass.render());
		loadGenomeSecondPassEntry.addCommand("touch $@");
		
		// add to make file
		make.addMakeEntry(loadGenomeSecondPassEntry);
		
		
		
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
			secondPassAlignEntry.setComment("Run second pass alignment for sample " + sample.getSampleId());
			secondPassAlignEntry.addDependency(genomeIndexSecondPassEntry);
			secondPassAlignEntry.addDependency(loadGenomeSecondPassEntry);
			
			secondPassAlignEntry.setTarget(dirBase +  "/" + sample.getSampleId()  + ".SECOND_PASS.OK");
			
			String sampleDir = dirBase + "/" + sample.getSampleId() + "/";
			
			ST SecondPassAlign = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndKeep --readFilesIn <files> --readFilesCommand <uncompress> --outFileNamePrefix <outdir>");
			SecondPassAlign.add("star", applicationOptions.getStar());
			SecondPassAlign.add("genomeDir", genomeOutDirPass2);
			SecondPassAlign.add("files", StringUtils.collectionToDelimitedString(sample.getFastqFiles(), " "));
			SecondPassAlign.add("uncompress", applicationOptions.getUncompressCommand());
			SecondPassAlign.add("outdir", sampleDir);
			
			
			
			
			secondPassAlignEntry.addCommand("mkdir -p " + sampleDir);
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
		removeGenomeEntry.addDependency(genomeIndexSecondPassEntry);
		removeGenomeEntry.addDependency(loadGenomeSecondPassEntry);
		
		removeGenomeEntry.addDependencies(pass2Dependencies);
		
		removeGenomeEntry.setTarget(baseDir + "/GENOME_REMOVED_SECOND_PASS.OK");
		
		ST removeGenomeSecondPass = new ST("<star> --genomeDir <genomeDir> --genomeLoad Remove");
		removeGenomeSecondPass.add("star", applicationOptions.getStar());
		removeGenomeSecondPass.add("genomeDir", genomeOutDirPass2);
		
		
		removeGenomeEntry.addCommand(removeGenomeSecondPass.render());
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
