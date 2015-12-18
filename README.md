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

# Download

```
wget https://github.com/christopher-gillies/rna-seq-pipeline/blob/master/release/rna-seq-pipeline-1.0.0.jar?raw=true
```
___

# Align fastq files with Star, find unique mapping reads, and mark duplicates
* Download STAR: https://github.com/alexdobin/STAR/archive/2.5.0b.zip
* Download Picard tools:  https://github.com/broadinstitute/picard/releases/download/1.141/picard-tools-1.141.zip
* Download Reference gene annotation: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz (http://feb2014.archive.ensembl.org/info/data/ftp/index.html)
* Download Reference genome: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
* The command below will generate a makefile in the specified output directory. The makefile contains the necessary commands for genome indexing, alignment, filtering, and duplicate marking
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
--readLength should be the lenght of the reads in your fastq files
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
* $OUT/SAMPLE_ID/final.bam in folder for each sample
###FORMAT 
```
[SAMPLE_ID][TAB][PATH_TO_BAM][TAB][LOG_PASS_1][TAB][LOG_PASS_2][SPLICE_JUNCTION_TABLE_PASS_1][TAB][SPLICE_JUNCTION_TABLE_PASS_2][TAB][DUPLICATE_METRIC_FILE]
#EXAMPLE
25969	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969_1//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25969//duplicate.output.metrics
25968	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//final.bam	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968_1//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//Log.final.out	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968_1//SJ.out.tab /net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//SJ.out.tab	/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015//25968//duplicate.output.metrics
```
___

#Quantify transcript and gene expression using flux capacitor
* Download flux capacitor
```
wget http://sammeth.net/artifactory/barna/barna/barna.capacitor/1.6.1/flux-capacitor-1.6.1.tgz
```

```
export FLUX_MEM=4G
export OUT=/net/assembly/cgillies/data/NEPTUNE/RNA-Seq/11_24_2015/
export PIPELINE=/net/wonderland/home/cgillies/programs/rna-seq-pipeline/target/rna-seq-pipeline-0.0.1-SNAPSHOT.jar 
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
--fluxCapacitor path to the flux capacitor executable file
--numberOfThreadsForFlux the number of threads for each flux command
-j for make command is the number fo flux capacitor instances to run in parallel
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

#Other
banner generated at http://patorjk.com/software/taag/
