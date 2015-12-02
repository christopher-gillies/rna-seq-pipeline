package org.kidneyomics.gtf;

public class VersionTrimmer {
	static final String pattern = "[.][0-9]+";
	public static String trim(String in) {
		return in.replaceAll(pattern, "");
	}
}
