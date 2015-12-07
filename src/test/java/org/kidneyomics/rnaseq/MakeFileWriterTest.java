package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

import ch.qos.logback.classic.Logger;

import static org.mockito.Mockito.*;

public class MakeFileWriterTest {

	
	@Test
	public void testStarAlign() throws Exception {
		
		//Setup
		File fastq1 = new File("/tmp/fastq1.fastq.gz");
		File fastq2 = new File("/tmp/fastq2.fastq.gz");
		File fastq3 = new File("/tmp/fastq3.fastq.gz");
		File fastq4 = new File("/tmp/fastq4.fastq.gz");
		
		StringBuilder sb = new StringBuilder();
		sb.append("SAMPLE1\t/tmp/fastq1.fastq.gz\t/tmp/fastq2.fastq.gz\n");
		sb.append("SAMPLE2\t/tmp/fastq3.fastq.gz\t/tmp/fastq4.fastq.gz\n");
		
		File fastqFileList = new File("/tmp/fastqFileList.txt");
		
		FileUtils.write(fastqFileList, sb.toString());
		
		List<File> files = new LinkedList<File>();
		
		files.add(fastq1);
		files.add(fastq2);
		files.add(fastq3);
		files.add(fastq4);
		files.add(fastqFileList);
		
		File reference = new File("/tmp/reference.fa");
		File star = new File("/tmp/star");
		File gtf = new File("/tmp/genecode.gtf");
		
		
		files.add(reference);
		files.add(star);
		files.add(gtf);
		
		//Create files
		//needed for validate method
		for(File f : files) {
			f.createNewFile();
		}		
		
		
		
		//Mock
		
		ApplicationOptions applicationOptionsMock = mock(ApplicationOptions.class);
		
		when(applicationOptionsMock.getFastqFiles()).thenReturn("/tmp/fastqFileList.txt");
		
		when(applicationOptionsMock.getReferenceSequence()).thenReturn("/tmp/reference.fa");
		
		when(applicationOptionsMock.getGtf()).thenReturn("/tmp/genecode.gtf");
		
		when(applicationOptionsMock.getStar()).thenReturn("/tmp/star");
		
		when(applicationOptionsMock.getOutputDirectory()).thenReturn("/tmp/test");
		
		when(applicationOptionsMock.getNumThreadsGenomeIndex()).thenReturn("11");
		
		when(applicationOptionsMock.getUncompressCommand()).thenReturn("zcat");
		
		when(applicationOptionsMock.getSjdbOverhang()).thenReturn(40);
		
		when(applicationOptionsMock.getNumThreadsAlign()).thenReturn("2");
		
		when(applicationOptionsMock.getJarLocation()).thenReturn("/data/RNA-Seq/11_18_2015/rna-seq-pipeline/target/rna-seq-pipeline-0.0.1-SNAPSHOT.jar");
		
		
		when(applicationOptionsMock.isNoSharedMemory()).thenReturn(false);
		
		LoggerService loggerService = new LoggerService();
		
		MakeFileWriter makeFileWriter = new MakeFileWriter(loggerService);
		makeFileWriter.applicationOptions = applicationOptionsMock;
		makeFileWriter.writeMakeFile(Mode.ALIGN);
		
		File out = new File("/tmp/test/Makefile");
		
		assertTrue(out.exists());
		
		String fileText = FileUtils.readFileToString(out);
		
		loggerService.getLogger(this).info("\n" + fileText);
		
		
		Resource r = new ClassPathResource("MakefileTestAlign");
		String fileTextExpect = FileUtils.readFileToString(r.getFile());
		
		assertEquals(fileTextExpect, fileText);
		
		//cleanup
		for(File f : files) {
			f.delete();
		}
	}
	
	
	@Test
	public void testStarAlignNoSharedMemory() throws Exception {
		
		//Setup
		File fastq1 = new File("/tmp/fastq1.fastq.gz");
		File fastq2 = new File("/tmp/fastq2.fastq.gz");
		File fastq3 = new File("/tmp/fastq3.fastq.gz");
		File fastq4 = new File("/tmp/fastq4.fastq.gz");
		
		StringBuilder sb = new StringBuilder();
		sb.append("SAMPLE1\t/tmp/fastq1.fastq.gz\t/tmp/fastq2.fastq.gz\n");
		sb.append("SAMPLE2\t/tmp/fastq3.fastq.gz\t/tmp/fastq4.fastq.gz\n");
		
		File fastqFileList = new File("/tmp/fastqFileList.txt");
		
		FileUtils.write(fastqFileList, sb.toString());
		
		List<File> files = new LinkedList<File>();
		
		files.add(fastq1);
		files.add(fastq2);
		files.add(fastq3);
		files.add(fastq4);
		files.add(fastqFileList);
		
		File reference = new File("/tmp/reference.fa");
		File star = new File("/tmp/star");
		File gtf = new File("/tmp/genecode.gtf");
		
		
		files.add(reference);
		files.add(star);
		files.add(gtf);
		
		//Create files
		//needed for validate method
		for(File f : files) {
			f.createNewFile();
		}		
		
		
		
		//Mock
		
		ApplicationOptions applicationOptionsMock = mock(ApplicationOptions.class);
		
		when(applicationOptionsMock.getFastqFiles()).thenReturn("/tmp/fastqFileList.txt");
		
		when(applicationOptionsMock.getReferenceSequence()).thenReturn("/tmp/reference.fa");
		
		when(applicationOptionsMock.getGtf()).thenReturn("/tmp/genecode.gtf");
		
		when(applicationOptionsMock.getStar()).thenReturn("/tmp/star");
		
		when(applicationOptionsMock.getOutputDirectory()).thenReturn("/tmp/test");
		
		when(applicationOptionsMock.getNumThreadsGenomeIndex()).thenReturn("11");
		
		when(applicationOptionsMock.getUncompressCommand()).thenReturn("zcat");
		
		when(applicationOptionsMock.getSjdbOverhang()).thenReturn(40);
		
		when(applicationOptionsMock.getNumThreadsAlign()).thenReturn("2");
		
		when(applicationOptionsMock.getJarLocation()).thenReturn("/data/RNA-Seq/11_18_2015/rna-seq-pipeline/target/rna-seq-pipeline-0.0.1-SNAPSHOT.jar");
		
		when(applicationOptionsMock.isNoSharedMemory()).thenReturn(true);
		
		LoggerService loggerService = new LoggerService();
		
		MakeFileWriter makeFileWriter = new MakeFileWriter(loggerService);
		makeFileWriter.applicationOptions = applicationOptionsMock;
		makeFileWriter.writeMakeFile(Mode.ALIGN);
		
		File out = new File("/tmp/test/Makefile");
		
		assertTrue(out.exists());
		
		String fileText = FileUtils.readFileToString(out);
		
		loggerService.getLogger(this).info("\n" + fileText);
		
		
		Resource r = new ClassPathResource("MakefileTestAlignNoShared");
		String fileTextExpect = FileUtils.readFileToString(r.getFile());
		
		assertEquals(fileTextExpect, fileText);
		
		//cleanup
		for(File f : files) {
			f.delete();
		}
	}
	
	@Test
	public void testWriteFluxCapacitorCommands() throws Exception {
		
		LoggerService loggerService = new LoggerService();

		
		//Setup
			File bam1 = new File("/tmp/tmp1.bam");
			File bam2 = new File("/tmp/tmp2.bam");
			
			StringBuilder sb = new StringBuilder();
			sb.append("SAMPLE1\t/tmp/tmp1.bam\n");
			sb.append("SAMPLE2\t/tmp/tmp2.bam\n");
			
			File bamlist = new File("/tmp/bamlist.txt");
			
			FileUtils.write(bamlist, sb.toString());
			
			List<File> files = new LinkedList<File>();
			
			files.add(bam1);
			files.add(bam2);
			files.add(bamlist);

			File flux = new File("/tmp/flux");
			File gtf = new File("/tmp/genecode.gtf");
			
			files.add(flux);
			files.add(gtf);
			
			//Create files
			//needed for validate method
			for(File f : files) {
				loggerService.getLogger(this).info("PATH:\t" + f.getAbsolutePath());
				f.createNewFile();
			}
			
			
			
			ApplicationOptions applicationOptionsMock = mock(ApplicationOptions.class);
			
			when(applicationOptionsMock.getBamList()).thenReturn(bamlist.getAbsolutePath());
			
			
			when(applicationOptionsMock.getGtf()).thenReturn(gtf.getAbsolutePath());
			
			when(applicationOptionsMock.getFluxCapacitor()).thenReturn(flux.getAbsolutePath());
			
			when(applicationOptionsMock.getOutputDirectory()).thenReturn("/tmp/test2/");
			
			when(applicationOptionsMock.getJarLocation()).thenReturn("/data/RNA-Seq/11_18_2015/rna-seq-pipeline/target/rna-seq-pipeline-0.0.1-SNAPSHOT.jar");
			
			when(applicationOptionsMock.getNumThreadsFlux()).thenReturn("10");
			
			
			
			MakeFileWriter makeFileWriter = new MakeFileWriter(loggerService);
			makeFileWriter.applicationOptions = applicationOptionsMock;
			makeFileWriter.writeMakeFile(Mode.FLUX_CAPACITOR);
			File out = new File("/tmp/test2/Makefile");
			
			assertTrue(out.exists());
			
			String fileText = FileUtils.readFileToString(out);
			
			loggerService.getLogger(this).info("\n" + fileText);
			
			File gtfList = new File("/tmp/test2/gtf.list.txt");
			
			files.add(gtfList);
			
			assertTrue(gtfList.exists());
			
			String gtfListText = FileUtils.readFileToString(gtfList);
			loggerService.getLogger(this).info(gtfListText);
			assertTrue(gtfListText.equals("SAMPLE1\t/tmp/test2//SAMPLE1.gtf\nSAMPLE2\t/tmp/test2//SAMPLE2.gtf\n"));
			
			
			Resource r = new ClassPathResource("MakefileTestFlux");
			String fileTextExpect = FileUtils.readFileToString(r.getFile());
			
			assertEquals(fileTextExpect, fileText);
			
			//cleanup
			for(File f : files) {
				f.delete();
			}
	}
	
}
