package org.kidneyomics.gtf;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.biojava.nbio.genome.parsers.gff.Feature;
import org.biojava.nbio.genome.parsers.gff.Location;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

public class FindOverlappingExonsBetweenGenesTest {

	Logger logger = LoggerFactory.getLogger(FindOverlappingExonsBetweenGenesTest.class);
	
	@Test
	public void testTwoGenes() throws IOException {
		logger.info("testTwoGenes");
		Resource r = new ClassPathResource("genecode.v19.annotation.head.10000.gtf.gz");
		
		File gtf = r.getFile();
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		LinkedList<Feature> gene2 = new LinkedList<>();
		
		
		try(GTFReader reader = GTFReader.getGTFByFileNoEnsemblVersion(gtf)) {
			//Only read exons
			reader.addFilter(new ExonFilter()).addFilter(new RemoveRetainedIntronFilter());
			for(Feature f : reader) {
				String geneId = f.getAttribute("gene_id");
				if(geneId.equals("ENSG00000157916")) {
					gene1.add(f);
				} else if(geneId.equals("ENSG00000157911")) {
					gene2.add(f);
				}
			}
		}
		
		GTFFeatureUtil.sortFeatures(gene1);
		GTFFeatureUtil.sortFeatures(gene2);
		
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertTrue(results.size() > 0);
		
		logger.info("Size of results " + results.size());
		
		gene1.removeAll(results);
		gene2.removeAll(results);
		
		List<Feature> gene1Merge = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(gene1);
		List<Feature> gene2Merge = FeatureMerger.mergeOverlappingFeaturesIgnoringStrand(gene2);
		
		LinkedList<Feature> noOverlap = new LinkedList<>();
		
		noOverlap.addAll(gene1Merge);
		noOverlap.addAll(gene2Merge);
		
		assertTrue(GTFFeatureUtil.hasNoOverlapIgnoreStrand(noOverlap));
		
	}
	
	
	@Test
	public void test1() {
		logger.info("test1");
		
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		List<Feature> results = service.findOverlappingExons(gene1, gene1);
		
		assertEquals(2,results.size());
	}
	
	
	@Test
	public void test2() {
		logger.info("test2");
		
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		Feature f3 = new Feature("chr2", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr2", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();
		gene2.add(f3);
		gene2.add(f4);
		
		Collections.sort(gene2, new FeatureComparator());
		
		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(0,results.size());
	}
	
	@Test
	public void test3() {
		logger.info("test3");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(199, 250, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr1", "a", "exon", Location.fromBio(400, 500, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();
		gene2.add(f3);
		gene2.add(f4);
		
		Collections.sort(gene2, new FeatureComparator());
		
		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(2,results.size());
	}
	
	@Test
	public void test4() {
		logger.info("test4");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(149, 250, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr1", "a", "exon", Location.fromBio(400, 500, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();
		gene2.add(f3);
		gene2.add(f4);
		
		Collections.sort(gene2, new FeatureComparator());
		
		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(3,results.size());
	}

	
	@Test
	public void test5() {
		logger.info("test5");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(149, 250, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr1", "a", "exon", Location.fromBio(45, 60, '+'), 0.0, 0, "");
		Feature f5 = new Feature("chr1", "a", "exon", Location.fromBio(37, 51, '+'), 0.0, 0, "");
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();
		gene2.add(f3);
		gene2.add(f4);
		gene2.add(f5);
		
		Collections.sort(gene2, new FeatureComparator());
		//for(Feature f : gene2) {
		//	System.err.println(GTFFeatureRenderer.render(f));
		//}
		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(5,results.size());
	}
	
	@Test
	public void test6() {
		logger.info("test6");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		Feature f6 = new Feature("chr1", "a", "exon", Location.fromBio(1, 37, '+'), 0.0, 0, "");
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		gene1.add(f6);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		Feature f3 = new Feature("chr1", "a", "exon", Location.fromBio(149, 250, '+'), 0.0, 0, "");
		Feature f4 = new Feature("chr1", "a", "exon", Location.fromBio(45, 60, '+'), 0.0, 0, "");
		Feature f5 = new Feature("chr1", "a", "exon", Location.fromBio(37, 51, '+'), 0.0, 0, "");
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();
		gene2.add(f3);
		gene2.add(f4);
		gene2.add(f5);
		
		Collections.sort(gene2, new FeatureComparator());
		//for(Feature f : gene2) {
		//	System.err.println(GTFFeatureRenderer.render(f));
		//}
		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(6,results.size());
	}
	
	@Test
	public void test7() {
		logger.info("test7");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		
		Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(100, 200, '+'), 0.0, 0, "");
		Feature f2 = new Feature("chr1", "a", "exon", Location.fromBio(50, 150, '+'), 0.0, 0, "");
		Feature f6 = new Feature("chr1", "a", "exon", Location.fromBio(1, 37, '+'), 0.0, 0, "");
		
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		gene1.add(f1);
		gene1.add(f2);
		gene1.add(f6);
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();

		
		Collections.sort(gene2, new FeatureComparator());
		//for(Feature f : gene2) {
		//	System.err.println(GTFFeatureRenderer.render(f));
		//}
		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(0,results.size());
	}
	
	@Test
	public void test8() {
		logger.info("test8");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		LinkedList<Feature> gene1 = new LinkedList<>();
		for(int i = 1000; i >= 1; i--) {
			Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(i+1, i + 100, '+'), 0.0, 0, "");
			gene1.add(f1);
		}
		
		Collections.sort(gene1, new FeatureComparator());
		
		
		LinkedList<Feature> gene2 = new LinkedList<>();
		for(int i = 1; i <= 1000; i++) {
			Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(i+1, i + 100, '+'), 0.0, 0, "");
			gene2.add(f1);
		}
		
		
		Collections.sort(gene2, new FeatureComparator());

		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(2000,results.size());
	}
	
	@Test
	public void test9() {
		logger.info("test9");
		FindOverlappingExonsBetweenGenes service = new FindOverlappingExonsBetweenGenes();
		
		/*
		 * 
> a = c()
> for(i in 1:1000) { if(i %% 3 == 0) { a = c(a,i) }  }
> a
  [1]   3   6   9  12  15  18  21  24  27  30  33  36  39  42  45  48  51  54  57  60  63  66  69  72  75  78  81  84  87  90  93  96  99 102 105 108 111 114 117 120 123 126 129 132 135 138
 [47] 141 144 147 150 153 156 159 162 165 168 171 174 177 180 183 186 189 192 195 198 201 204 207 210 213 216 219 222 225 228 231 234 237 240 243 246 249 252 255 258 261 264 267 270 273 276
 [93] 279 282 285 288 291 294 297 300 303 306 309 312 315 318 321 324 327 330 333 336 339 342 345 348 351 354 357 360 363 366 369 372 375 378 381 384 387 390 393 396 399 402 405 408 411 414
[139] 417 420 423 426 429 432 435 438 441 444 447 450 453 456 459 462 465 468 471 474 477 480 483 486 489 492 495 498 501 504 507 510 513 516 519 522 525 528 531 534 537 540 543 546 549 552
[185] 555 558 561 564 567 570 573 576 579 582 585 588 591 594 597 600 603 606 609 612 615 618 621 624 627 630 633 636 639 642 645 648 651 654 657 660 663 666 669 672 675 678 681 684 687 690
[231] 693 696 699 702 705 708 711 714 717 720 723 726 729 732 735 738 741 744 747 750 753 756 759 762 765 768 771 774 777 780 783 786 789 792 795 798 801 804 807 810 813 816 819 822 825 828
[277] 831 834 837 840 843 846 849 852 855 858 861 864 867 870 873 876 879 882 885 888 891 894 897 900 903 906 909 912 915 918 921 924 927 930 933 936 939 942 945 948 951 954 957 960 963 966
[323] 969 972 975 978 981 984 987 990 993 996 999
> length(a)
[1] 333
		 */
		LinkedList<Feature> gene1 = new LinkedList<>();
		for(int i = 1000; i >= 1; i--) {
			if(i % 3 == 0) {
				Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(i, i, '+'), 0.0, 0, "");
				gene1.add(f1);
			}
		}
		
		Collections.sort(gene1, new FeatureComparator());
		
/*
> b=c()
> for(i in 1:1000) { if(i %% 6 == 0) { b = c(b,i) }  }
> b
  [1]   6  12  18  24  30  36  42  48  54  60  66  72  78  84  90  96 102 108 114 120 126 132 138 144 150 156 162 168 174 180 186 192 198 204 210 216 222 228 234 240 246 252 258 264 270 276
 [47] 282 288 294 300 306 312 318 324 330 336 342 348 354 360 366 372 378 384 390 396 402 408 414 420 426 432 438 444 450 456 462 468 474 480 486 492 498 504 510 516 522 528 534 540 546 552
 [93] 558 564 570 576 582 588 594 600 606 612 618 624 630 636 642 648 654 660 666 672 678 684 690 696 702 708 714 720 726 732 738 744 750 756 762 768 774 780 786 792 798 804 810 816 822 828
[139] 834 840 846 852 858 864 870 876 882 888 894 900 906 912 918 924 930 936 942 948 954 960 966 972 978 984 990 996
> length(b)
[1] 166
 */
 
		LinkedList<Feature> gene2 = new LinkedList<>();
		for(int i = 1; i <= 1000; i++) {
			if(i % 6 == 0) {
				Feature f1 = new Feature("chr1", "a", "exon", Location.fromBio(i, i, '+'), 0.0, 0, "");
				gene2.add(f1);
			}
		}
		
		
		Collections.sort(gene2, new FeatureComparator());

		
		List<Feature> results = service.findOverlappingExons(gene1, gene2);
		
		assertEquals(2 * 166,results.size());
	}
}
