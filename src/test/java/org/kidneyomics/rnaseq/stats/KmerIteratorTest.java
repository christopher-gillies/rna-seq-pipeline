package org.kidneyomics.rnaseq.stats;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.junit.Test;

public class KmerIteratorTest {

	@Test
	public void test() {
		String testData = "1234567890";
		int k = 3;
		KmerIterator iter = new KmerIterator(k, testData);
		//length = 10 - 3 + 1
		
		List<String> kmers = new ArrayList<String>();
		while(iter.hasNext()) {
			String kmer = iter.next();
			System.err.println(kmer);
			kmers.add(kmer);
		}
		
		assertEquals("123",kmers.get(0));
		assertEquals("234",kmers.get(1));
		assertEquals("345",kmers.get(2));
		assertEquals("456",kmers.get(3));
		assertEquals("567",kmers.get(4));
		assertEquals("678",kmers.get(5));
		assertEquals("789",kmers.get(6));
		assertEquals("890",kmers.get(7));
		
		assertEquals(8,kmers.size());
		
		assertEquals(null,iter.next());
		
	}
	
	@Test
	public void test2() {
		String testData = "12345678901";
		int k = 4;
		KmerIterator iter = new KmerIterator(k, testData);
		//length = 11 - 4 + 1
		
		List<String> kmers = new ArrayList<String>();
		while(iter.hasNext()) {
			String kmer = iter.next();
			System.err.println(kmer);
			kmers.add(kmer);
		}
		
		assertEquals("1234",kmers.get(0));
		assertEquals("2345",kmers.get(1));
		assertEquals("3456",kmers.get(2));
		assertEquals("4567",kmers.get(3));
		assertEquals("5678",kmers.get(4));
		assertEquals("6789",kmers.get(5));
		assertEquals("7890",kmers.get(6));
		assertEquals("8901",kmers.get(7));
		
		assertEquals(8,kmers.size());
		
		assertEquals(null,iter.next());
		
	}
	
	@Test
	public void test3() {
		String testData = "12345678901";
		int k = 2;
		KmerIterator iter = new KmerIterator(k, testData);
		//length = 11 - 2 + 1 = 10
		
		List<String> kmers = new ArrayList<String>();
		while(iter.hasNext()) {
			String kmer = iter.next();
			System.err.println(kmer);
			kmers.add(kmer);
		}
		
		assertEquals("12",kmers.get(0));
		assertEquals("23",kmers.get(1));
		assertEquals("34",kmers.get(2));
		assertEquals("45",kmers.get(3));
		assertEquals("56",kmers.get(4));
		assertEquals("67",kmers.get(5));
		assertEquals("78",kmers.get(6));
		assertEquals("89",kmers.get(7));
		assertEquals("90",kmers.get(8));
		assertEquals("01",kmers.get(9));
		
		assertEquals(10,kmers.size());
		
		assertEquals(null,iter.next());
		
	}
	
	@Test
	public void test4() {
		String testData = "12345678901";
		int k = 1;
		KmerIterator iter = new KmerIterator(k, testData);
		//length = 11 - 1 + 1 = 11
		
		List<String> kmers = new ArrayList<String>();
		while(iter.hasNext()) {
			String kmer = iter.next();
			System.err.println(kmer);
			kmers.add(kmer);
		}
		
		assertEquals("1",kmers.get(0));
		assertEquals("2",kmers.get(1));
		assertEquals("3",kmers.get(2));
		assertEquals("4",kmers.get(3));
		assertEquals("5",kmers.get(4));
		assertEquals("6",kmers.get(5));
		assertEquals("7",kmers.get(6));
		assertEquals("8",kmers.get(7));
		assertEquals("9",kmers.get(8));
		assertEquals("0",kmers.get(9));
		assertEquals("1",kmers.get(10));
		
		assertEquals(11,kmers.size());
		
		assertEquals(null,iter.next());
		
	}

	@Test
	public void test6() {
		String testData = "12345678901";
		int k = 11;
		KmerIterator iter = new KmerIterator(k, testData);
		//length = 11 - 11 + 1 = 11
		
		List<String> kmers = new ArrayList<String>();
		while(iter.hasNext()) {
			String kmer = iter.next();
			System.err.println(kmer);
			kmers.add(kmer);
		}
		
		assertEquals("12345678901",kmers.get(0));
		
		assertEquals(1,kmers.size());
		
		assertEquals(null,iter.next());
		
	}
}
