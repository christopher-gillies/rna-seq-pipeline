package org.kidneyomics.util;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Chr2Int {

	private static Pattern pattern = Pattern.compile("[0-9]+");
	
	public static int convert(String in) {
		in = in.replaceFirst("[Cc][hH][rR]", "");
		Matcher matcher = pattern.matcher(in);
		if(matcher.matches()) {
			return Integer.parseInt(in);
		} else {
			if(in.equalsIgnoreCase("X")) {
				return 23;
			} else if(in.equalsIgnoreCase("Y")) {
				return 24;
			} else if( in.equalsIgnoreCase("M")) {
				return 25;
			} else if( in.equalsIgnoreCase("MT")) {
				return 25;
			} else {
				return 26;
			}
		}
		
		
	}
}
