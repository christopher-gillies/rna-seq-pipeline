package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.junit.Test;

public class VersionTrimmerTest {

	@Test
	public void test() {
		String res = VersionTrimmer.trim("ENTS000001111.1");
		assertEquals("ENTS000001111",res);
	}

}
