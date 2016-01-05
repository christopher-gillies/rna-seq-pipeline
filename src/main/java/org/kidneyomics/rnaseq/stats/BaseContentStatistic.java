package org.kidneyomics.rnaseq.stats;

import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.SAMRecord;

class BaseContentStatistic extends AbstractReadPairStatistic {

	private long A;
	private long T;
	private long G;
	private long C;
	private long N;
	
	private static final byte A_BYTE_VAL = 'A';
	private static final byte T_BYTE_VAL = 'T';
	private static final byte G_BYTE_VAL = 'G';
	private static final byte C_BYTE_VAL = 'C';
	private static final byte N_BYTE_VAL = 'N';
	
	// "A\tT\tG\tC\tN\tGC\tGT\tAT\tAC";
	private static final String[] fields = { "A", "T", "G", "C", "N", "GC", "GT", "AT", "AC" };
	
	@Override
	protected void addRecord(SAMRecord record) {
		final byte[] read = record.getReadBases();
		if(record.getReadNegativeStrandFlag()) {
			//A=T , G≡C
			for(int i = 0; i < read.length; i++) {
				switch(read[i]) {
				case A_BYTE_VAL:
					T++;
					break;
				case T_BYTE_VAL:
					A++;
					break;
				case G_BYTE_VAL:
					C++;
					break;
				case C_BYTE_VAL:
					G++;
					break;
				case N_BYTE_VAL:
					N++;
					break;
					default:
						throw new RuntimeException("Unsupported base in read: " + record.getReadString());
				}
			}
		} else {
			for(int i = 0; i < read.length; i++) {
				switch(read[i]) {
				case A_BYTE_VAL:
					A++;
					break;
				case T_BYTE_VAL:
					T++;
					break;
				case G_BYTE_VAL:
					G++;
					break;
				case C_BYTE_VAL:
					C++;
					break;
				case N_BYTE_VAL:
					N++;
					break;
					default:
						throw new RuntimeException("Unsupported base in read: " + record.getReadString());
				}
			}
		}
	}

	@Override
	public List<Double> getStatistic() {
		double totalBases = A + T + G + C + N;
		LinkedList<Double> result = new LinkedList<>();
		/*
		System.err.println("A: " + A);
		System.err.println("T: " + T);
		System.err.println("G: " + G);
		System.err.println("C: " + C);
		System.err.println("N: " + N);
		System.err.println("Total: " + totalBases);
		*/
		//A
		result.add(  A / totalBases );
		//T
		result.add(  T / totalBases );
		//G
		result.add(  G / totalBases );
		//C
		result.add(  C / totalBases );
		//N
		result.add(  N / totalBases );
		
		//G C
		result.add(  (G + C) / totalBases );
		//G T
		result.add(  (G + T) / totalBases );
		//A T
		result.add(  (A + T) / totalBases );
		//A C
		result.add(  (A + C) / totalBases );
		
		return result;
	}

	@Override
	protected String[] getFields() {
		return fields;
	}

}
