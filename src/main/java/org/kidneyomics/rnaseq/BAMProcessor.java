package org.kidneyomics.rnaseq;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BAMProcessor implements AutoCloseable {

	private SamReader reader;
	private SAMRecordIterator iterator;
	private Queue<SAMRecord> queue;
	
	private BAMProcessor(File in) {
		reader = SamReaderFactory.makeDefault().open(in);
		iterator = reader.iterator();
		queue = new LinkedList<SAMRecord>();
	}
	
	public static BAMProcessor getBAMProcessor(File in) {
		return new BAMProcessor(in);
	}
	
	/**
	 * This is needed when bam is sorted in coordinate order instead of query order
	 * @return the next pair of mapped reads, if only one read is mapped then return that one read
	 * return null if no mapped reads are left.
	 */
	public SAMRecordPair getNextReadPair() {
		SAMRecord mate1;
		//Get first mapped read
		while( (mate1 = getFirstMate()) != null && mate1.getReadUnmappedFlag() == true); 
		if(mate1 == null) {
			return null;
		}
		
		//At this point mate1 should not be null;
		if(mate1.getMateUnmappedFlag()) {
			SAMRecordPair srp = new SAMRecordPair();
			srp.setMate1(mate1);
			return srp;
		} 
		
		//find mate
		SAMRecord mate2 = getSecondMate(mate1);
		if(mate2 == null) {
			throw new IllegalStateException("Second mate not found for:\n" + mate1.toString());
		} else {
			SAMRecordPair srp = new SAMRecordPair();
			srp.setMate1(mate1);
			srp.setMate2(mate2);
			return srp;
		}
		
	}
	
	SAMRecord getFirstMate() {
		//check queue, and then check the record iterator
		SAMRecord mate1 = null;
		if(queue.peek() != null) {
			mate1 = queue.poll();
		} else {
			if(iterator.hasNext()) {
				mate1 = iterator.next();
			}
		}
		return mate1;
		
	}
	
	SAMRecord getSecondMate(SAMRecord mate1) {
		
		SAMRecord mate2 = null;
		Iterator<SAMRecord> queueIter = queue.iterator();
		while(mate2 == null && queueIter.hasNext()) {
			SAMRecord tmp = queueIter.next();
			if(areMates(mate1,tmp)) {
				mate2 = tmp;
				queueIter.remove();
			}
		}
		
		if(mate2 != null) {
			return mate2;
		} 
		
		while(mate2 == null && iterator.hasNext()) {
			SAMRecord tmp = iterator.next();
			if(areMates(mate1,tmp)) {
				mate2 = tmp;
			} else {
				//if the record is unmapped
				if(!tmp.getReadUnmappedFlag()) {
					queue.add(tmp);
				}
			}
		}
		
		return mate2;
		
		
	}
	
	/*
	 * Maybe relax this to only ensure that the names match?
	 */
	public boolean areMates(SAMRecord mate1, SAMRecord mate2) {
		if(mate1.getMateAlignmentStart() == mate2.getAlignmentStart()
				&& mate2.getMateAlignmentStart() == mate1.getAlignmentStart() 
				&& mate1.getReadName().equals(mate2.getReadName())) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public void close() throws Exception {
		iterator.close();	
	}

	
}
