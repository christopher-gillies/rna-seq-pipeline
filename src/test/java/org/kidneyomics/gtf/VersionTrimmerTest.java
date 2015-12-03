package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import org.junit.Test;

public class VersionTrimmerTest {

	@Test
	public void test1() {
		String res = VersionTrimmer.trim("ENST000001111.1");
		assertEquals("ENST000001111",res);
	}

	@Test
	public void test2() {
		String res = VersionTrimmer.trim("ENSG000001111.1");
		assertEquals("ENSG000001111",res);
	}
	
	
	
	@Test
	public void test3() {
		String res = VersionTrimmer.trim("ENSG000001111.1\tENSG000001111.1\t53.1\tCCCC");
		assertEquals("ENSG000001111\tENSG000001111\t53.1\tCCCC",res);
	}
}
