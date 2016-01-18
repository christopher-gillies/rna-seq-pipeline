```
 ___  _ _  ___       ___  ___  ___  
| . \| \ || . | ___ / __>| __>| . | 
|   /|   ||   ||___|\__ \| _> | | | 
|_\_\|_\_||_|_|     <___/|___>`___\ 
                                    
 ___  _  ___  ___  _    _  _ _  ___ 
| . \| || . \| __>| |  | || \ || __>
|  _/| ||  _/| _> | |_ | ||   || _> 
|_|  |_||_|  |___>|___||_||_\_||___>
```

[![Build Status](https://travis-ci.org/christopher-gillies/rna-seq-pipeline.svg?branch=master)](https://travis-ci.org/christopher-gillies/rna-seq-pipeline?branch=master)

* This program generates makefiles to process RNA-seq data. The methodology is similar to the GTEx methodology. There are a few differences. For example this program uses STAR aligner instead of Tophat.
* This program starts with a list of fastq files and creates a series of BAM files using STAR aligner
* The BAM files can then be quantified using flux capacitor or a read counting program. Both these programs output a series of gene expression matrices in terms of counts and RPKM that can be analyzed further.
* To align samples you will need to download STAR aligner, picard tools, a reference gene annotation, and a reference genome. Currently, Ensembl and GENCODE are supported.
* To quantify transcript expression you will need to install flux capacitor
* To quantify exon and gene expression you just need the bam list produced from the alignment step
* This program is implemented in Java using Spring Boot, Apache commons, HTSJDK, ANTLR String template and biojava 4.1.
* Please make sure you have java 1.7 or later installed for this program to work

___
# Download

```
wget -O rna-seq-pipeline-1.0.0.jar https://github.com/christopher-gillies/rna-seq-pipeline/blob/master/release/rna-seq-pipeline-1.0.0.jar?raw=true
```
___
___
# Align fastq files with Star, find unique mapping reads, and mark duplicates
* Download STAR: https://github.com/alexdobin/STAR/archive/2.5.0b.zip
* Download Picard tools:  https://github.com/broadinstitute/picard/releases/download/1.141/picard-tools-1.141.zip
* Download Reference gene annotation: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz (http://feb2014.archive.ensembl.org/info/data/ftp/index.html)
* Download Reference genome: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
* The command below will generate a makefile in the specified output directory. The makefile contains the necessary commands for genome indexing, alignment, filtering, and duplicate marking
* This pipeline follows the 2-pass approach mentioned in the GATK best practices: https://www.broadinstitute.org/gatk/guide/article?id=3891
```
export PIPELINE=/path/to/rna-seq-pipeline/rna-seq-pipeline-1.0.0.jar
export STAR=/data/RNA-Seq/11_18_2015/STAR/bin/Linux_x86_64_static/STAR
export OUT=/data/RNA-Seq/11_18_2015/
export GTF=$OUT/Homo_sapiens.GRCh37.75.gtf
export REF=$OUT/Homo_sapiens.GRCh37.75.dna.primary_assembly.no.extra.contigs.fa
export PICARD=/home/cgillies/FluidigmReferences/programs/picard-tools-1.135/picard.jar
export FASTQ_LIST=/data/RNA-Seq/11_18_2015/fastq.files.txt

#Unzip gtf and assembly file
zcat $OUT/Homo_sapiens.GRCh37.75.gtf.gz > $OUT/Homo_sapiens.GRCh37.75.gtf
zcat $OUT/Homo_sapiens.GRCh37.75.dna.primary_assembly.no.extra.contigs.fa.gz > $OUT/Homo_sapiens.GRCh37.75.dna.primary_assembly.no.extra.contigs.fa


java -jar $PIPELINE --fastqList $FASTQ_LIST --gtf $GTF --numberOfThreadsForGenomeIndex 16 --outputDir $OUT --referenceSequence $REF --starAligner $STAR --readLength 78 --numberOfThreadsForAlign 2 --picard $PICARD

cd $OUT

nohup make -j 2 1> nohup.out 2> nohup.err &
```

## Parameter notes
```
--readLength should be the length of the reads in your fastq files
--numberOfThreadsForGenomeIndex parameter should be as many threads as you can spare because genome indexing is time consuming. This is done one time and before all samples are aligned.
--numberOfThreadsForAlign is the number of threads to use for the align command
-j with make command specifies the number of samples to align at the same time
```

```
--fastqList
#FORMAT
[SAMPLE_ID]TAB[PATH_TO_FASTQ_MATE_1]TAB[PATH_TO_FASTQ_MATE_2]

#EXAMPLE
25848	/path/R1.clipped.fastq.gz	/path/R2.clipped.fastq.gz
25848	/path/R1.clipped.fastq.gz	/path/R2.clipped.fastq.gz
25848	/path/R1.clipped.fastq.gz	/path/R2.clipped.fastq.gz
```
##Files in output directory
* $OUT/bam.list.txt -- a list of bam files
* $OUT/SAMPLE_ID/final.bam -- in folder for each sample
* $OUT/STAR_RUN_STATS.txt -- a table containing all the statistics from the final STAR log for each sample
* $OUT/MERGED_DUP_STATS.txt -- a table summarizing the results from the mar duplicate scripts
* $OUT/MERGED_BAM_STATS.txt -- a table summarizing basic statistics across samples such as GC content, mean phred score, mean Q30 bases per gene
###Format of bam.list.txt
```
[SAMPLE_ID][TAB][PATH_TO_BAM][TAB][LOG_PASS_1][TAB][LOG_PASS_2][SPLICE_JUNCTION_TABLE_PASS_1][TAB][SPLICE_JUNCTION_TABLE_PASS_2][TAB][DUPLICATE_METRIC_FILE]
#EXAMPLE
25969	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//duplicate.output.metrics
25968	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968_1//SJ.out.tab /net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//duplicate.output.metrics
```
___
___

#Quantify transcript and gene expression using flux capacitor
* this program will generate a makefile to quantify sample expression using flux capacitor. It will also merge the expression values into a series gene expression matrices. 
* download flux capacitor
```
wget http://sammeth.net/artifactory/barna/barna/barna.capacitor/1.6.1/flux-capacitor-1.6.1.tgz
```

```
export FLUX_MEM=4G
export OUT=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/
export PIPELINE=/path/to/rna-seq-pipeline/rna-seq-pipeline-1.0.0.jar
export FLUX_CAPACITOR=~/programs/flux-capacitor-1.6.1/bin/flux-capacitor
export GTF=$OUT/Homo_sapiens.GRCh37.75.gtf
export BAMLIST=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/bam.list.txt 
export OUTDIR=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/FLUX/
java -jar $PIPELINE --fluxCapacitor $FLUX_CAPACITOR --bamList $BAMLIST --gtf $GTF --outputDir $OUTDIR --numberOfThreadsForFlux 5
cd $OUT
nohup make -j2 1> nohup.out 2> nohup.err &
```

##Parameter information
```
--bamList is the bam list from the alignment command
--fluxCapacitor path to the flux capacitor executable file
--numberOfThreadsForFlux the number of threads for each flux command
-j for make command is the number of flux capacitor instances to run in parallel
```
##Files in output directory
* sampleid.gtf -- separate gtf file containing reads and rpkm per transcript
* gene.count.expression.txt -- a matrix containing the sum of reads across all transcripts per sample
* gene.rpkm.expression.txt -- a matrix containing the sum of RPKM across all transcripts per sample
* transcript.count.expression.txt -- a matrix containing the transcript read counts per sample
* transcript.rpkm.expression.txt -- a matrix containing the transcript read rpkm per sample
* transcript.ratios.from.count.expression.txt -- a matrix containing the transcript ratio using transcript count data per sample. The ratio for a transcript in a subject is the transcript read count / sum of transcripts read counts for the gene of the transcript of interest.
* transcript.ratios.from.rpkm.expression.txt -- a matrix containing the transcript ratio using transcript rpkm data per sample. The ratio for a transcript in a subject is the transcript rpkm / sum of all transcript rpkms for the gene of the transcript of interest.
* gtf.list.txt -- a tab separated file containing the sample id and the path to the gtf file for the sample of interest

___
___

#Count reads in exons and genes using GTEx-like methodology
* currently this program is written for unstranded, paired-end RNA-seq data
* these read counts follow closely to the methodology mentioned in the GTEx paper http://www.ncbi.nlm.nih.gov/pubmed/25954001
* exons labeled 'retained_intron' are excluded
* exons that overlap across genes are excluded
* remaining exons that overlap are merged on a per gene basis
* the program counts read pairs not unpaired reads
* reads have to be uniquely mapped
* reads have to be properly paired
* alignment edit distance has to be <= 6
* read pairs or split reads that map to multiple genes are excluded
* reads pairs overlapping multiple exons from the same gene are fractionally counted. For example if a has 78 mapped bases, where 30 mapped to one exon, 20 to a second exon and 38 to a third exon, then the exons would be counted 30/78, 20/78, and 38/78 respectively.
* One difference between the GTEx methodology and this program is that this program fractionally counts reads that overlap introns, where GTEx ignores reads where this happens. GTEx requires reads to be 100% within exon boundaries. For example, if a read has 78 mapped bases, and 32 overlap an exon, then this read will be counted by this program as 32/78.
* Gene-level counts are calculated as the sum of exons for the gene. The length of the gene is equal to the sum of the length of the exons included in the analysis. So exons that overlap between genes are not included in the calculation of the gene length, because no reads are counted over such regions.
* rpkm values are normalized by the number of uniquely mapped reads because these are the only reads used in quantification.
```
export OUT=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/
export OUT_DIR=$OUT/EXON_COUNTS/
export BAM_LIST=$OUT/bam.list.txt
export GTF=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/Homo_sapiens.GRCh37.75.gtf
export PIPELINE=/path/to/rna-seq-pipeline/rna-seq-pipeline-1.0.0.jar
java -jar $PIPELINE --countReadsInAllSamples --bamList $BAM_LIST --gtf $GTF --outputDir $OUT_DIR
cd $OUT_DIR
nohup make -j10 1> nohup.out 2> nohup.err &
```

##Parameter information
```
--bamList is the bam list from the alignment command
```
##Output files
* gene.exon.counts.rpkm.sample_id.gtf -- gtf file containing the exon and gene rpkms and read counts
* gene.exon.counts.rpkm.sample_id.gtf.stats -- statistics file for each gtf. This file contains the total number of read pairs, the number and fraction of ambiguously mapped reads, the number and fraction of reads filtered out, the number and fraction of unmapped reads, the number and fraction of partially mapped reads, and finally the number and fraction of mapped reads
* exon.counts.txt -- a count expression matrix for each exon and sample
* exon.rpkm.txt	 -- a rpkm expression matrix for each exon and sample
* gene.counts.txt -- a count expression matrix for each gene and sample
* gene.rpkm.txt -- a rpkm expression matrix for each gene and sample
* merged.stats.txt -- the statistics files merged for each sample

___
___

#Find uniquely mapped reads in a STAR-aligned BAM
* to use this function, just download the release jar file and specify the input and output bam file names as described below
```
export PIPELINE=/path/to/rna-seq-pipeline/rna-seq-pipeline-1.0.0.jar
export FILE_IN=/path/to/bam
export FILE_OUT=/path/to/newbam
java -jar $PIPELINE --findUniqueMappedReads --fileIn $FILE_IN --fileOut $FILE_OUT
```
___
___

#Relabel expression matrix
* this function is useful if you need to change the sample ids from the expression matrix to new values
* you need to specify the expression matrix that you want to relabel (--expressionMatrix)
* you need to specify the new path of the expression matrix after relabeling (--fileOut)
* lastly, you need to specify a two column tab seprated file (--fileIn), where first column is the old id and the second column is the new id. Each line is a separate mapping.
* This command works by constructing a dictionary mapping between the old ids and the new ids. If you only specify a single id mapping, then this command will only relabel the single sample and all other ids will remain the same.

```
export OUT=/OUT_DIRECTORY/
export IN_FILE=/PATH/ID_MAP_FOR_EXPRESSION.txt
export EXPRESSION_MATRIX=/PATH/gene.rpkm.txt
export OUT_FILE=$OUT/EXON_COUNTS//reheader.gene.rpkm.txt
export PIPELINE=/PATH/rna-seq-pipeline-1.0.0.jar
java -jar $PIPELINE --mapExpressionIds --fileIn $IN_FILE --expressionMatrix $EXPRESSION_MATRIX --fileOut $OUT_FILE
```
##Here is the --fileIn format for this command
```
OLD_ID1 NEW_ID1
OLD_ID2 NEW_ID2
OLD_ID3 NEW_ID3
OLD_ID4 NEW_ID4
```
___
___

#Other
banner generated at http://patorjk.com/software/taag/
