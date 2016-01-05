package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import org.junit.Test;

public class StreamingMeanVarianceCalculatorTest {

	@Test
	public void test() {
		StreamingMeanAndVarianceCalculator calc = new StreamingMeanAndVarianceCalculator();
		
		calc.add(17.0);
		calc.add(19.0);
		calc.add(24.0);
		
		assertEquals(20,calc.getMean(),0.00001);
		assertEquals(13,calc.getVariance(),0.00001);
		assertEquals(3.605551,calc.getSd(),0.00001);
	}
	
	
	@Test
	public void test2() {
		StreamingMeanAndVarianceCalculator calc = new StreamingMeanAndVarianceCalculator();
		
		for(int i = 1; i <= 100; i++) {
			calc.add(i);
		}
		
		assertEquals(50.5,calc.getMean(),0.00001);
		assertEquals(841.6667,calc.getVariance(),0.001);
		assertEquals(29.01149,calc.getSd(),0.001);
	}

}
