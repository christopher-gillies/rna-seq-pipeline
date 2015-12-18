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


# Align with Star
* Download STAR: https://github.com/alexdobin/STAR/archive/2.5.0b.zip
* Download Picard tools:  https://github.com/broadinstitute/picard/releases/download/1.141/picard-tools-1.141.zip
* Download Reference gene annotation: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz (http://feb2014.archive.ensembl.org/info/data/ftp/index.html)
* Download Reference genome: ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

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
###FORMAT 
```
[SAMPLE_ID][TAB][PATH_TO_BAM][SPLICE_JUNCTION_TABLE_PASS_1][TAB][SPLICE_JUNCTION_TABLE_PASS_2][TAB][DUPLICATE_METRIC_FILE]
```







#Other
banner generated at http://patorjk.com/software/taag/
