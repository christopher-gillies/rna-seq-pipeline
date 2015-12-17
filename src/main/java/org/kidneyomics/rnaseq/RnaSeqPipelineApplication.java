package org.kidneyomics.rnaseq;




import java.io.File;

import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.ApplicationContext;
import org.springframework.context.annotation.Import;

@SpringBootApplication
@Import(RnaSeqPipelineConfiguration.class)
public class RnaSeqPipelineApplication {

	

    public static void main(String[] args) throws Exception {
    	    	
    	
        //ApplicationContext context = SpringApplication.run(RnaSeqPipelineApplication.class, args);
        
    	SpringApplication springApplication = new SpringApplication(new Object[] { RnaSeqPipelineApplication.class });
    	springApplication.setLogStartupInfo(false);
    	ApplicationContext context = springApplication.run(args);

        

    	Logger logger = LoggerFactory.getLogger(RnaSeqPipelineApplication.class);

    	
        
        ApplicationOptions applicationOptions = context.getBean(ApplicationOptions.class);
		
       
        
        
        logger.info("Application location: " + applicationOptions.getJarLocation());
        
        
        try {
	        Mode mode = applicationOptions.validateOptions();
			
	        switch(mode) {
	        
	        case ALIGN: {
	        	logger.info("Creating makefile for alignment");
	        	MakeFileWriter makeFileWriter = context.getBean(MakeFileWriter.class);
	        	makeFileWriter.writeMakeFile(mode);
	        	break;
	        }
	        case FIND_UNIQUE_MAPPED_READS: {
	        	logger.info("Removing multimapped reads");
	        	UniqueMappingFilter filter = context.getBean(UniqueMappingFilter.class);
	        	File in = new File(applicationOptions.getFileIn());
	        	File out = new File(applicationOptions.getFileOut());
	        	filter.filter(in,out);
	        	break;
	        }
	        case FLUX_CAPACITOR: {
	        	logger.info("Creating makefile for transcript deconvolution");
	        	MakeFileWriter makeFileWriter = context.getBean(MakeFileWriter.class);
	        	makeFileWriter.writeMakeFile(mode);
	        	break;
	        }
	        case GENE_EXPRESSION_MATRIX: {
	        	logger.info("Creating gene expression matrix");
	        	FluxMerge fluxMerge = context.getBean(FluxMerge.class);
	        	fluxMerge.writeGeneMatrix();
	        	break;
	        }
	        case TRANSCRIPT_EXPRESSION_MATRIX: {
	        	logger.info("Creating transcript expression matrix");
	        	FluxMerge fluxMerge = context.getBean(FluxMerge.class);
	        	fluxMerge.writeTranscriptMatrix();
	        	break;
	        }
	        case TRANSCRIPT_RATIO_MATRIX:
	        	logger.info("Creating transcript ratio matrix");
	        	FluxMerge fluxMerge = context.getBean(FluxMerge.class);
	        	fluxMerge.writeTranscriptRatioMatrix();
	        	break;
	        case COUNT_READS_IN_EXONS: {
	        	logger.info("Counting read pairs in exons");
	        	ExonQuantifier exonQuantifier = context.getBean(ExonQuantifier.class);
	        	exonQuantifier.quantify();
	        	break;
	        }
	        case MERGE_EXON_COUNTS: {
	        	logger.info("Merging exon gtf count files");
	        	ExonMerge exm = context.getBean(ExonMerge.class);
	        	exm.writeOutMatrices();
	        	break;
	        }
	        case COUNT_READS_ALL_SAMPLES: {
	        	logger.info("Creating makefile to count all reads across samples");
	        	MakeFileWriter makeFileWriter = context.getBean(MakeFileWriter.class);
	        	makeFileWriter.writeMakeFile(mode);
	        	break;
	        }
	        case ERROR: {
	        	
	        }
	        }
        } catch(Exception e) {
        	logger.error(e.getMessage());
        	throw e;
        }
		
		
    }
    
    
}
