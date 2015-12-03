package org.kidneyomics.gtf;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VersionTrimmer {
	static //static final String pattern = "[.][0-9]+";
	Pattern pattern = Pattern.compile("(ENS[TGE][0-9]+)[.][0-9]+");
	public static String trim(String in) {
		Matcher m = pattern.matcher(in);
		if(m.find()) {
			return m.replaceAll("$1");
		} else {
			return in;
		}
		
	}
}
