package org.kidneyomics.rnaseq;




import java.io.File;

import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.ApplicationContext;

@SpringBootApplication
public class RnaSeqPipelineApplication {

	

    public static void main(String[] args) throws Exception {
    	    	
    	
        //ApplicationContext context = SpringApplication.run(RnaSeqPipelineApplication.class, args);
        
    	SpringApplication springApplication = new SpringApplication(new Object[] { RnaSeqPipelineApplication.class });
    	springApplication.setLogStartupInfo(false);
    	ApplicationContext context = springApplication.run(args);

        

    	Logger logger = LoggerFactory.getLogger(RnaSeqPipelineApplication.class);

    	
        
        ApplicationOptions applicationOptions = context.getBean(ApplicationOptions.class);
		
       
        
        
        logger.info("Application location: " + applicationOptions.getJarLocation());
        
        Mode mode = applicationOptions.validateOptions();
		
        switch(mode) {
        
        case ALIGN: {
        	logger.info("Creating makefile for alignment");
        	MakeFileWriter makeFileWriter = context.getBean(MakeFileWriter.class);
        	makeFileWriter.writeMakeFile(mode);
        	break;
        }
        case FIND_UNIQUE_MAPPED_READS:
        	logger.info("Removing multimapped reads");
        	UniqueMappingFilter filter = context.getBean(UniqueMappingFilter.class);
        	File in = new File(applicationOptions.getFileIn());
        	File out = new File(applicationOptions.getFileOut());
        	filter.filter(in,out);
        	
        	break;
        case ERROR: {
        	
        }
        }
		
		
    }
    
    
}
