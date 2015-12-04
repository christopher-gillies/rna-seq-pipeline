package org.kidneyomics.gtf;

import java.util.Comparator;

import org.biojava.nbio.genome.parsers.gff.Location;

public class LocationComparator implements Comparator<Location> {

	@Override
	public int compare(Location o1, Location o2) {
		return Integer.compare(o1.bioStart(), o2.bioStart());
	}

}
