package org.kidneyomics.rnaseq;




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
			ApplicationCommand command = null;
	        
	        switch(mode) {
	        
	        case COUNT_READS_ALL_SAMPLES:
	        case FLUX_CAPACITOR:
	        case ALIGN: {
	        	logger.info("Creating makefile for alignment");
	        	command = context.getBean(MakeFileWriter.class);
	        	break;
	        }
	        case FIND_UNIQUE_MAPPED_READS: {
	        	logger.info("Removing multimapped reads");
	        	command = context.getBean(UniqueMappingFilter.class);
	        	break;
	        }
	        case TRANSCRIPT_EXPRESSION_MATRIX:
	        case TRANSCRIPT_RATIO_MATRIX:
	        case GENE_EXPRESSION_MATRIX: {
	        	logger.info("Running flux merge");
	        	command = context.getBean(FluxMerge.class);
	        	break;
	        }
	        case COUNT_READS_IN_EXONS: {
	        	logger.info("Counting read pairs in exons");
	        	command = context.getBean(ExonQuantifier.class);
	        	break;
	        }
	        case MERGE_EXON_COUNTS: {
	        	logger.info("Merging exon gtf count files");
	        	command = context.getBean(ExonMerge.class);
	        	break;
	        }
	        case MERGE_EXON_COUNTS_STATS: {
	        	logger.info("Merging statistics for all exon counts");
	        	command = context.getBean(ReadCountStatMerger.class);
	        	break;
	        }
	        case MERGE_STAR_LOGS: {
	        	logger.info("Merging statistics for all star logs");
	        	command = context.getBean(STARLogMerger.class);
	        	break;
	        }
	        case MAP_EXPRESSION_IDS:
	        	logger.info("Remapping expression ids");
	        	command = context.getBean(MapIdentifiersForExpressionMatrix.class);
	        	break;
	        case HELP:
	        	//logger.info("Showing help");
	        	command = context.getBean(DoNothingApplicationCommand.class);
	        case ERROR: {
	        	
	        }
	        }
	        
	        if(command != null) {
	        	command.doWork();
	        } else {
	        	throw new Exception("No command selected");
	        }
        } catch(Exception e) {
        	logger.error(e.getMessage());
        	throw e;
        }
		
		
    }
    
    
}
