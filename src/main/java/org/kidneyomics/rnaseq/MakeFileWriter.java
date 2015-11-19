package org.kidneyomics.rnaseq;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

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
		
		Map<String,SampleData> sampleMap = new HashMap<String,SampleData>();
		
		
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
			
			ST firstPassAlign = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndKeep --readFilesIn <files> --readFilesCommand <uncompress> --outFileNamePrefix <outdir> --outSJfilterCountUniqueMin 4 2 2 2 --outSJfilterCountTotalMin 4 2 2 2 --runThreadN <n>");
			firstPassAlign.add("star", applicationOptions.getStar());
			firstPassAlign.add("genomeDir", genomeOutDir);
			firstPassAlign.add("files", StringUtils.collectionToDelimitedString(sample.getFastqFiles(), " "));
			firstPassAlign.add("uncompress", applicationOptions.getUncompressCommand());
			firstPassAlign.add("outdir", sampleDir);
			firstPassAlign.add("n", applicationOptions.getNumThreadsAlign());
			
			String sjdbFile = sampleDir + "/SJ.out.tab";
			String samFile = sampleDir + "/Aligned.out.sam";
			String finalLog = sampleDir + "/Log.final.out";
			sjdbFiles.add(sjdbFile);
			
			
			firstPassAlignEntry.addCommand("mkdir -p " + sampleDir);
			firstPassAlignEntry.addCommand(firstPassAlign.render());
			firstPassAlignEntry.addCommand("touch $@");
			
			// add to make file
			make.addMakeEntry(firstPassAlignEntry);
			
			//add to dependencies
			pass1Dependencies.add(firstPassAlignEntry);
			
			
			/*
			 * Create sample entry in map
			 */
			SampleData sampleData = new SampleData();
			sampleData.id = sample.getSampleId();
			sampleData.pass1Dir = sampleDir;
			sampleData.samPass1 = samFile;
			sampleData.finalLogPass1 = finalLog;
			sampleData.sjdb1 = sjdbFile;
			
			sampleMap.put(sampleData.id, sampleData);
			
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
			
			ST SecondPassAlign = new ST("<star> --genomeDir <genomeDir> --genomeLoad LoadAndKeep --readFilesIn <files> --readFilesCommand <uncompress> --outFileNamePrefix <outdir> --runThreadN <n>");
			SecondPassAlign.add("star", applicationOptions.getStar());
			SecondPassAlign.add("genomeDir", genomeOutDirPass2);
			SecondPassAlign.add("files", StringUtils.collectionToDelimitedString(sample.getFastqFiles(), " "));
			SecondPassAlign.add("uncompress", applicationOptions.getUncompressCommand());
			SecondPassAlign.add("outdir", sampleDir);
			SecondPassAlign.add("n", applicationOptions.getNumThreadsAlign());
			
			
			
			secondPassAlignEntry.addCommand("mkdir -p " + sampleDir);
			secondPassAlignEntry.addCommand(SecondPassAlign.render());
			secondPassAlignEntry.addCommand("touch $@");
			
			String sjdbFile = sampleDir + "/SJ.out.tab";
			String samFile = sampleDir + "/Aligned.out.sam";
			String finalLog = sampleDir + "/Log.final.out";
			// add to make file
			make.addMakeEntry(secondPassAlignEntry);
			
			//add to dependencies
			pass2Dependencies.add(secondPassAlignEntry);
			
			SampleData sampleData = sampleMap.get(sample.getSampleId());
			if(sampleData == null) {
				throw new NullPointerException(sample.getSampleId() + " is missing...");
			} else {
				sampleData.samPass2 = samFile;
				sampleData.pass2Dir = sampleDir;
				sampleData.sjdb2 = sjdbFile;
				sampleData.finalLogPass2 = finalLog;
				sampleData.polishDependency = secondPassAlignEntry;
			}
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
		 *  AddOrReplaceReadGroups and MarkDuplicates
		 *  "Polish Bam"
		 * 
		 */
		
		for(SampleData sampleData : sampleMap.values()) {
			MakeEntry polishBamEntry = new MakeEntry();
			sampleData.cleanUpDependency = polishBamEntry;
			polishBamEntry.setComment("Add the readgroups to the bam and mark duplicates for " + sampleData.id);
			polishBamEntry.setTarget(dirBase + "/SAMPLE_" + sampleData.id + "_BAM_POLISH.OK");
			polishBamEntry.addDependency(sampleData.polishDependency);
			
			sampleData.sortedBam = sampleData.pass2Dir + "/sorted.bam";
			
			
			ST addOrReplaceReadGroups = new ST("java -jar <picard> AddOrReplaceReadGroups I=<input> O=<output> SO=coordinate RGID=<rgid> RGLB=<library> RGPL=ILLUMINA RGPU=ILLUMINA RGSM=<sample>");
			addOrReplaceReadGroups.add("picard", applicationOptions.getPicard())
			.add("input", sampleData.samPass2)
			.add("output", sampleData.sortedBam)
			.add("rgid", sampleData.id)
			.add("library", sampleData.id)
			.add("sample", sampleData.id);
			
			polishBamEntry.addCommand(addOrReplaceReadGroups.render());
			
			sampleData.finalBam = sampleData.pass2Dir + "/final.bam";
			sampleData.dupMetrics = sampleData.pass2Dir + "/duplicate.output.metrics";
			ST markDuplicatesEntry = new ST("java -jar <picard> MarkDuplicates I=<input> O=<output>  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=<metrics>");
			markDuplicatesEntry.add("picard", applicationOptions.getPicard())
			.add("input", sampleData.sortedBam)
			.add("output", sampleData.finalBam)
			.add("metrics", sampleData.dupMetrics);
			
			polishBamEntry.addCommand(markDuplicatesEntry.render());
			
			
			polishBamEntry.addCommand("touch $@");
			
			make.addMakeEntry(polishBamEntry);
		}
		
		
		
		
		/*
		 * 
		 * 
		 * Remove unncecessary files
		 * 
		 * 
		 * 
		 */
		
		for(SampleData sampleData : sampleMap.values()) {
			MakeEntry cleanUpEntry = new MakeEntry();

			cleanUpEntry.setComment("Remove unnessary files for " + sampleData.id);
			cleanUpEntry.setTarget(dirBase + "/SAMPLE_" + sampleData.id + "_CLEAN_UP.OK");
			cleanUpEntry.addDependency(sampleData.cleanUpDependency);
			
			cleanUpEntry.addCommand("rm " + sampleData.samPass1);
			cleanUpEntry.addCommand("rm " + sampleData.samPass2);
			cleanUpEntry.addCommand("rm " + sampleData.sortedBam);
			cleanUpEntry.addCommand("touch $@");
			
			make.addMakeEntry(cleanUpEntry);
		}
		
		
		/*
		 * 
		 * 
		 * Write output data
		 * 
		 * 
		 * 
		 */
		
		StringBuilder sb = new StringBuilder();
		for(SampleData sampleData : sampleMap.values()) {
			sb.append(sampleData.toLine());
			sb.append("\n");
		}
		
		FileUtils.write(new File(dirBase + "bam.list.txt"), sb.toString());
		
		
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
	
	
	private class SampleData {
		String id;
		String pass1Dir;
		String pass2Dir;
		String finalLogPass1;
		String finalLogPass2;
		String samPass1;
		String samPass2;
		String sortedBam;
		String finalBam;
		String sjdb1;
		String sjdb2;
		String dupMetrics;
		MakeEntry polishDependency;
		MakeEntry cleanUpDependency;
		
		public String toLine() {
			return StringUtils.arrayToDelimitedString(new String[] { id, finalBam, finalLogPass1, finalLogPass2, sjdb1, sjdb2, dupMetrics  }, "\t");
		}
	}
	
}
