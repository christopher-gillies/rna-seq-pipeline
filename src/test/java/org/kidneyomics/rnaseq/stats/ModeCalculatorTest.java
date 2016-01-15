package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import org.junit.Test;

public class ModeCalculatorTest {

	@Test
	public void test() {
		ModeCalculator calc = new ModeCalculator();
		
		for(int i = 0; i <= 100; i++) {
			for(int j = 0; j < i + 1; j++) {
				calc.add(i);
			}
		}
		
		
		assertEquals(100,calc.getMode());
		
		for(int j = 0; j < 1000; j++) {
			calc.add(50);
		}
		
		assertEquals(50,calc.getMode());
	}

}
