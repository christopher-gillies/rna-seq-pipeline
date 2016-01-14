package org.kidneyomics.rnaseq;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Test;

public class DelimitedFileEntryTest {

	@Test
	public void test1() {
		DelimitedFileEntry entry = DelimitedFileEntry.getDelimitedFileEntry("Sample1", new String[] { "col1","col2", "col3" },new String[] { "val1", "val2", "val3"} );
		
		assertEquals("val1",entry.getValue("col1"));
		assertEquals("val2",entry.getValue("col2"));
		assertEquals("val3",entry.getValue("col3"));
		
		assertEquals("NA",entry.getValue("x"));
		
		//assertEquals("col1\tcol2\tcol3",DelimitedFileEntry.headerToString(cols, delimiter))
		
	}
	
	
	@Test
	public void test2() throws IOException {
		DelimitedFileEntry entry = DelimitedFileEntry.getDelimitedFileEntry("Sample1", new String[] { "col1","col2", "col3" },new String[] { "val11", "val12", "val13"} );
		DelimitedFileEntry entry2 = DelimitedFileEntry.getDelimitedFileEntry("Sample2", new String[] { "col1","col2", "col3","col4" },new String[] { "val21", "val22", "val23","val24"} );
		
		List<DelimitedFileEntry> list = new LinkedList<>();
		
		list.add(entry);
		list.add(entry2);
		
		List<String> header = DelimitedFileEntry.getAllKeys(list);
		
		assertEquals("SAMPLE_ID\tcol1\tcol2\tcol3\tcol4",DelimitedFileEntry.headerToString(header, "\t"));
		
		assertEquals("Sample1\tval11\tval12\tval13\tNA",entry.toString(header, "\t"));
		assertEquals("Sample2\tval21\tval22\tval23\tval24",entry2.toString(header, "\t"));
		
		String tmpDir = FileUtils.getTempDirectoryPath();
		
		File out = new File(tmpDir + "/test.txt");
		
		DelimitedFileEntry.writeToFile(list, out.getAbsolutePath(), "\t");
		
		assertTrue(out.exists());
		
		List<String> lines = FileUtils.readLines(out);
		assertEquals(3,lines.size());
		
		assertEquals("SAMPLE_ID\tcol1\tcol2\tcol3\tcol4",lines.get(0));
		assertEquals("Sample1\tval11\tval12\tval13\tNA",lines.get(1));
		assertEquals("Sample2\tval21\tval22\tval23\tval24",lines.get(2));
		
		out.delete();
	}

	@Test
	public void test3() {
		DelimitedFileEntry entry = DelimitedFileEntry.getDelimitedFileEntryFromDelimitedLine("S1", "col1\tcol2\tcol3\tcol4", "val21\tval22\tval23\tval24", "\t");
		
		assertEquals("val21",entry.getValue("col1"));
		assertEquals("val22",entry.getValue("col2"));
		assertEquals("val23",entry.getValue("col3"));
		assertEquals("val24",entry.getValue("col4"));
		
	}
}
