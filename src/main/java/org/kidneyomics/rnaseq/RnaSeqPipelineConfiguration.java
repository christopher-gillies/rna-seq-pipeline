package org.kidneyomics.rnaseq;

import org.kidneyomics.gtf.FindOverlappingFeatures;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Scope;

@Configuration
public class RnaSeqPipelineConfiguration {
	
	@Bean()
	@Scope("prototype")
	public FindOverlappingFeatures findOverlappingFeatures() {
		return new FindOverlappingFeatures();
	}
	
	
	@Bean()
	public QuantificationFactory quantificationFactory() {
		return new QuantificationFactory();
	}
}
