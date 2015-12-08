package org.kidneyomics.gtf;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.springframework.core.convert.converter.Converter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

import java.util.LinkedList;
import java.util.List;

public class SAMRecordToFeatureConverter implements Converter<SAMRecord, List<Feature>> {

	@Override
	public List<Feature> convert(SAMRecord in) {
		
		List<Feature> list = new LinkedList<Feature>();
		
		Cigar cigar = in.getCigar();
		List<CigarElement> elements = cigar.getCigarElements();
		elements.get(0).getOperator();
		elements.get(0).getLength();
		
		
		
		return list;
	}

	
	
}
