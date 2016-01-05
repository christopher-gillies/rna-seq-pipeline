package org.kidneyomics.rnaseq.stats;

class StreamingMeanAndVarianceCalculator {

	//http://www.johndcook.com/blog/standard_deviation/
	
	private long k = 0;
	
	private double mean = 0;
	private double S = 0;
	
	
	double getMean() {
		return mean;
	}
	
	double getVariance() {
		if(k < 2) {
			return 0;
		} else {
			return S / ( (double) k - 1);
		}
	}
	
	double getSd() {
		if(k < 2) {
			return 0;
		} else {
			return Math.sqrt(getVariance());
		}
	}
	
	void add(double val) {
		double oldMean = mean;
		if(k == 0) {
			mean = val;
		} else {
			mean = oldMean + (val - oldMean) / ( (double) k + 1 );
		}
		
		if(k < 1) {
			S = 0;
		} else if(k >= 1) {
			S = S + (val - oldMean)*(val - mean);
		}
		
		k++;
	}
	
}
