package org.kidneyomics.rnaseq;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.stereotype.Component;

@Component
public class LoggerService {
	
	
	public Logger getLogger(Object instance) {
		return LoggerFactory.getLogger(instance.getClass());
	}
	
}
