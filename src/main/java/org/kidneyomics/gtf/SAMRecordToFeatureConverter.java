package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.convert.converter.Converter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.LinkedList;
import java.util.List;

public class SAMRecordToFeatureConverter implements Converter<SAMRecord, List<Feature>> {

	Logger logger = LoggerFactory.getLogger(SAMRecordToFeatureConverter.class);
	
	@Override
	public List<Feature> convert(SAMRecord in) {
		
		
		/*
		 * example
		 * 
		 * 10M5N7M
		 * 
		 * start = 100
		 * 
		 * MMMMMMMMMMNNNNNMMMMMMM
		 * 
		 * F1 = start = 100, end = 109,
		 * F2 = start = 115, end = 121
		 * last element is inclusive
		 * 
		 * start = 100
		 * end = 100 + 10 - 1 = 109
		 * 
		 * newStart = 109 + 1 + 5 = 115
		 */
		
		List<Feature> list = new LinkedList<Feature>();
		
		Cigar cigar = in.getCigar();
		
		List<CigarElement> elements = cigar.getCigarElements();
		
		boolean isNegativeStrand = in.getReadNegativeStrandFlag();
		
		
		int currentStart = in.getAlignmentStart();
		int currentLength = 0;
		CigarElement previous = null;
		for(final CigarElement element : elements) {
			//logger.info(element.getLength() + "");
            switch (element.getOperator()) {
            case M:
            case EQ:
            case X:
            	currentLength += element.getLength();
                break;
            case D:
            case N:
            	if(previous.getOperator().equals(CigarOperator.D) || previous.getOperator().equals(CigarOperator.N)) {
            		//Just keep skipping if the previous element was a D or N
            		currentStart += element.getLength();
            	} else {
	            	int currentEnd = currentStart + currentLength - 1;
	            	
	            	list.add( buildFeature(in,currentStart,currentEnd,isNegativeStrand)  );
	            	
	            	
	            	currentStart = currentEnd + element.getLength() + 1;
	            	currentLength = 0;
            	}
            	break;
            }
            previous = element;
		}
		
		//add last feature
		int currentEnd = currentStart + currentLength - 1;
		list.add( buildFeature(in,currentStart,currentEnd,isNegativeStrand)  );
		
		return list;
	}
	
	public int getNumberOfMappedBases(SAMRecord in) {
		Cigar cigar = in.getCigar();
		List<CigarElement> elements = cigar.getCigarElements();
		int mappedBases = 0;
		for(final CigarElement element : elements) {
			//logger.info(element.getLength() + "");
            switch (element.getOperator()) {
            case M:
            	mappedBases += element.getLength();
                break;
            default:
            }
		}
		return mappedBases;
	}
	
	private Feature buildFeature(SAMRecord in, int start, int end, boolean isNegativeStrand) {
		Feature f = new Feature(in.getReferenceName(),"","",Location.fromBio(start, end, isNegativeStrand ? '-' : '+'),Double.MIN_VALUE,Integer.MAX_VALUE,"");
		return f;
	}

	
	
}
