package org.kidneyomics.rnaseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component()
public class STARLogMerger implements ApplicationCommand {

	//STAR
	
	/*
cgillies@assembly:/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/25986$ cat Log.final.out
                                 Started job on |	Nov 24 13:11:41
                             Started mapping on |	Nov 24 13:12:15
                                    Finished on |	Nov 24 13:18:34
       Mapping speed, Million of reads per hour |	303.09

                          Number of input reads |	31908652
                      Average input read length |	78
                                    UNIQUE READS:
                   Uniquely mapped reads number |	27319614
                        Uniquely mapped reads % |	85.62%
                          Average mapped length |	77.67
                       Number of splices: Total |	7769525
            Number of splices: Annotated (sjdb) |	7757722
                       Number of splices: GT/AG |	7616623
                       Number of splices: GC/AG |	124419
                       Number of splices: AT/AC |	19037
               Number of splices: Non-canonical |	9446
                      Mismatch rate per base, % |	0.22%
                         Deletion rate per base |	0.00%
                        Deletion average length |	1.43
                        Insertion rate per base |	0.00%
                       Insertion average length |	1.17
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |	2644721
             % of reads mapped to multiple loci |	8.29%
        Number of reads mapped to too many loci |	18984
             % of reads mapped to too many loci |	0.06%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |	0.00%
                 % of reads unmapped: too short |	5.96%
                     % of reads unmapped: other |	0.07%
	 */
	
	
	//DUPLICATE STATS
	
	/*
cgillies@assembly:/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/25986$ cat duplicate.output.metrics
## htsjdk.samtools.metrics.StringHeader
# picard.sam.markduplicates.MarkDuplicates INPUT=[/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25986//sorted.bam] OUTPUT=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/25986/final.bam METRICS_FILE=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/25986/duplicate.output.metrics VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true    MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 SORTING_COLLECTION_SIZE_RATIO=0.25 PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_NAME=MarkDuplicates REMOVE_DUPLICATES=false ASSUME_SORTED=false DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES READ_NAME_REGEX=[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).* OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 VERBOSITY=INFO QUIET=false COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json
## htsjdk.samtools.metrics.StringHeader
# Started on: Tue Nov 24 20:11:39 EST 2015

## METRICS CLASS	picard.sam.DuplicationMetrics
LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
25986	0	27319614	0	0	15168986	5531061	0.555242	16648213

## HISTOGRAM	java.lang.Double
BIN	VALUE
1.0	1.104634
2.0	1.318698
3.0	1.360181
4.0	1.36822
5.0	1.369778
6.0	1.37008
7.0	1.370138
8.0	1.37015
9.0	1.370152
10.0	1.370152
11.0	1.370152
12.0	1.370152
13.0	1.370152
14.0	1.370152
15.0	1.370152
16.0	1.370152
17.0	1.370152
18.0	1.370152
19.0	1.370152
20.0	1.370152
21.0	1.370152
22.0	1.370152
23.0	1.370152
24.0	1.370152
25.0	1.370152
26.0	1.370152
27.0	1.370152
28.0	1.370152
29.0	1.370152
30.0	1.370152
31.0	1.370152
32.0	1.370152
33.0	1.370152
34.0	1.370152
35.0	1.370152
36.0	1.370152
37.0	1.370152
38.0	1.370152
39.0	1.370152
40.0	1.370152
41.0	1.370152
42.0	1.370152
43.0	1.370152
44.0	1.370152
45.0	1.370152
46.0	1.370152
47.0	1.370152
48.0	1.370152
49.0	1.370152
50.0	1.370152
51.0	1.370152
52.0	1.370152
53.0	1.370152
54.0	1.370152
55.0	1.370152
56.0	1.370152
57.0	1.370152
58.0	1.370152
59.0	1.370152
60.0	1.370152
61.0	1.370152
62.0	1.370152
63.0	1.370152
64.0	1.370152
65.0	1.370152
66.0	1.370152
67.0	1.370152
68.0	1.370152
69.0	1.370152
70.0	1.370152
71.0	1.370152
72.0	1.370152
73.0	1.370152
74.0	1.370152
75.0	1.370152
76.0	1.370152
77.0	1.370152
78.0	1.370152
79.0	1.370152
80.0	1.370152
81.0	1.370152
82.0	1.370152
83.0	1.370152
84.0	1.370152
85.0	1.370152
86.0	1.370152
87.0	1.370152
88.0	1.370152
89.0	1.370152
90.0	1.370152
91.0	1.370152
92.0	1.370152
93.0	1.370152
94.0	1.370152
95.0	1.370152
96.0	1.370152
97.0	1.370152
98.0	1.370152
99.0	1.370152
100.0	1.370152
	 */
	
	
	/*
	 * BAM LIST
	 */
	
	/*
25969	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//duplicate.output.metrics
25968	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//duplicate.output.metrics
25963	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25963//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25963_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25963//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25963_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25963//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25963//duplicate.output.metrics
25962	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25962//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25962_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25962//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25962_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25962//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25962//duplicate.output.metrics
25961	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25961//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25961_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25961//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25961_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25961//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25961//duplicate.output.metrics
25960	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25960//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25960_1//Log.final.out	/net/assembly
	 */
	
	
	private ApplicationOptions applicationOptions;
	
	private Logger logger;
	
	@Autowired
	STARLogMerger(LoggerService loggerService, ApplicationOptions applicationOptions, SampleGTFReader sampleGtfReader, QuantificationFactory quantificationFactory) {
		this.logger = loggerService.getLogger(this);
		this.applicationOptions = applicationOptions;
	}
	
	public void readBamList() throws Exception {
		
		File bamList = new File(applicationOptions.getBamList());
		List<String> lines = FileUtils.readLines(bamList);
		logger.info("Reading bam list");
		List<BAMInfo> infos = BAMInfo.getBAMInfoFromLines(lines);
			
		List<STARLog> logs = new LinkedList<>();
		
		logger.info("Parsing star logs");
		for(BAMInfo info : infos) {
			String filePath = info.getLog2();
			logger.info("Reading log for " + info.getSampleId());
			List<String> logLines = FileUtils.readLines(new File(filePath));
			
			//process lines
			STARLog log = STARLog.createLogFromLines(info.getSampleId(),logLines);
			logs.add(log);
		
		}
		
		//get all keys
		Set<String> allKeys = new HashSet<>();
		for(STARLog log : logs) {
			allKeys.addAll(log.getKeys());
		}
		
		//sort the keys in lexicographic ordering
		List<String> header = new LinkedList<>();
		header.addAll(allKeys);
		Collections.sort(header);
		
		//Create writer
		String outfile = applicationOptions.getFileOut();
		Path p = Paths.get(outfile);
		BufferedWriter bf = Files.newBufferedWriter(p,Charset.defaultCharset());
		
		
		logger.info("Writing to " + outfile);
		//Write header
		bf.write("SampleId");
		bf.write("\t");
		bf.write(StringUtils.join(header, '\t'));
		bf.write("\n");
		
		
		//Write lines
		for(STARLog log : logs) {
			bf.write(log.getId());
			bf.write("\t");
			Iterator<String> iter = header.iterator();
			while(iter.hasNext()) {
				String key = iter.next();
				String value = log.getValue(key);
				bf.write(value);
				if(iter.hasNext()) {
					bf.write("\t");
				}
			}
			bf.write("\n");
		}
		
		//close
		bf.close();
		
		logger.info("finished");
	}

	@Override
	public void doWork() throws Exception {
		readBamList();
	}
}
