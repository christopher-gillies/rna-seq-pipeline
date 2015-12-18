package org.kidneyomics.rnaseq;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Queue;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BAMProcessorImplDict implements BAMProcessor {

	
	private HashMap<String,SAMRecordPair> readBuffer = new HashMap<>();
	
	
	private SamReader reader;
	private SAMRecordIterator iterator;
	private Queue<SAMRecord> queue;
	private long readCount = 0;
	private Logger logger = LoggerFactory.getLogger(BAMProcessorImplQueue.class);
	private int logSkipSize = 100000;
	
	
	public BAMProcessorImplDict withLogSkipSize(int logSkipSize) {
		this.logSkipSize = logSkipSize;
		return this;
	}
	
	private BAMProcessorImplDict(File in) {
		reader = SamReaderFactory.makeDefault().open(in);
		iterator = reader.iterator();
		queue = new LinkedList<SAMRecord>();
	}
	
	public static BAMProcessorImplDict getBAMProcessor(File in) {
		return new BAMProcessorImplDict(in);
	}
	
	@Override
	public SAMRecordPair getNextReadPair() {
		//insert first read into dictionary by queryname
		//insert second read into dictionary
		//check if the dictionary length for that entry has both pairs
		//if it is does return the read pair
		//otherwise continue reading
		//this way just return pairs as they are completed
		//should be MUCH faster
		//make sure to delete the entry after returning so that we dont have a memory leak
		
		if(iterator.hasNext()) {
			while(iterator.hasNext()) {
				SAMRecord record = iterator.next();
				countRead(record);
				//skip if the read is unmapped, not properly paired or mate is unmapped
				if(record.getReadUnmappedFlag() == true || record.getProperPairFlag() == false || record.getMateUnmappedFlag() == true ) {
					continue;
				}
				
				String query = record.getReadName();
				
				//check if read mate has been read already
				if(readBuffer.containsKey(query)) {
					//if it has then return the pair
					SAMRecordPair pair = readBuffer.get(query);
					pair.addPair(record);
					if(pair.bothPairsAligned() && pair.isValidPair()) {
						//prevent memory leak by deleting keys that are no longer needed
						readBuffer.remove(query);
						return pair;
					} else {
						throw new RuntimeException(query + " is not properly mated");
					}
				} else {
					//otherwise create an entry and store it by its query name
					SAMRecordPair pair = new SAMRecordPair();
					pair.addPair(record);
					readBuffer.put(query, pair);
				}
			}
		} else {
			if(readBuffer.size() > 0) {
				for(String key : readBuffer.keySet()) {
					logger.info("No mate for for " + key);
				}
				throw new RuntimeException("No mates found for some reads please make sure all reads are properly paired");
			}
		}
		
		return null;
	}
	

	private void countRead(final SAMRecord record) {
		if(record != null) {
			readCount++;
			if(readCount % logSkipSize == 0) {
				logger.info("Single reads processed: " + readCount);
				logger.info("Last read: " + record.getReferenceName() + ":" + record.getAlignmentStart() + "-" + record.getAlignmentEnd());
			}
		}
	}
	
	@Override
	public void close() throws Exception {
		iterator.close();	
	}
}
