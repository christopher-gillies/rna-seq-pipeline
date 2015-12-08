package org.kidneyomics.rnaseq;

import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ExonQuantifier {

	
	@Autowired
	ApplicationOptions applicationOptions;
	
	Logger logger;
	
	@Autowired
	ExonQuantifier(LoggerService loggerService) {
		this.logger = loggerService.getLogger(this);
	}
	
	
	
}
