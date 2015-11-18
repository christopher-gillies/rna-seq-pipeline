package org.kidneyomics.rnaseq;



import org.kidneyomics.rnaseq.ApplicationOptions.Mode;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import org.springframework.context.ApplicationContext;

@SpringBootApplication
public class RnaSeqPipelineApplication {

	

    public static void main(String[] args) throws Exception {
    	    	
        ApplicationContext context = SpringApplication.run(RnaSeqPipelineApplication.class, args);
        
        OptionProcessor op = context.getBean(ApplicationOptionProcessor.class);
        op.processInputs(args);
        
        ApplicationOptions applicationOptions = context.getBean(ApplicationOptions.class);
		Mode mode = applicationOptions.validateOptions();
		
		/*
		 * Perform more logic
		 */
		
		MakeFileWriter makeFileWriter = context.getBean(MakeFileWriter.class);
		
		makeFileWriter.writeMakeFile(mode);
		
    }
    
    
}
