# NGS3 DNA-seq: Pre-processing WGS data (Upstream analysis)

## [TABLE OF CONTENTS](#table-of-contents)
- [NGS3 DNA-seq: Pre-processing WGS data (Upstream analysis)](#ngs3-dna-seq-pre-processing-wgs-data-upstream-analysis)
  - [TABLE OF CONTENTS](#table-of-contents)
- [1. Introduction](#1-introduction)
- [2. Raw Data Processing](#2-raw-data-processing)
  - [Introduce to FastQ files](#introduce-to-fastq-files)
    - [File format](#file-format)
  - [Sequencing quality control (FastQC)](#sequencing-quality-control-fastqc)
    - [🌟 Hands-on: Check FastQC](#-hands-on-check-fastqc)
  - [Read trimming and filtering (Trimmomatic)](#read-trimming-and-filtering-trimmomatic)
    - [🌟 Hands-on: Trimming raw data](#-hands-on-trimming-raw-data)
- [3. Alignment: Mapping reads to reference genome](#3-alignment-mapping-reads-to-reference-genome)
  - [Introduce to SAM/BAM file](#introduce-to-sambam-file)
    - [Align with BWA mem](#align-with-bwa-mem)
    - [Convert SAM to BAM](#convert-sam-to-bam)
- [4. Mapped reads post-processing](#4-mapped-reads-post-processing)
  - [🌟 \[Hands-on\]: Use GATK to analyze data](#-hands-on-use-gatk-to-analyze-data)
  - [Sorting and marking duplicate, inxdexing the BAM file](#sorting-and-marking-duplicate-inxdexing-the-bam-file)
    - [SortSam](#sortsam)
    - [MarkDuplicates](#markduplicates)
    - [Collect Insert Size](#collect-insert-size)
- [5. Alignment Data: Quality Control](#5-alignment-data-quality-control)
  - [Coverage Analysis](#coverage-analysis)
    - [Case 1: Keep duplicated reads](#case-1-keep-duplicated-reads)
    - [Case 2: Remove duplicated reads](#case-2-remove-duplicated-reads)
  - [A coverage plot in R](#a-coverage-plot-in-r)
    - [Case 1](#case-1)
    - [Case 2](#case-2)

~

# 1. Introduction
In DNA sequencing, upstream analysis refers to the bioinformatic analysis that is performed on the raw sequencing data to process, filter, and prepare the data for downstream analysis.

# 2. Raw Data Processing

Create shortcuts to work dir

```bash
p_raw='/path/to/raw_file/'
p_trim='path/to/trimmed_file/'
p_align='path/to/aligned/file/'
p_ref='/path/to/references/'
```

## Introduce to FastQ files
https://en.wikipedia.org/wiki/FASTQ_format
### File format
![fastq format](https://www.drive5.com/usearch/manual9.2/fastq_fig.jpg)



## Sequencing quality control (FastQC)
### 🌟 Hands-on: Check FastQC

Check the quality of raw fastq file using FastQC:

```bash
cd $p_raw

# Basic syntax
fastqc <file> -o <path/to/output>

# Do for all file
for i in $(ls *.fastq.gz); do fastqc $i; done
```
Hands-on: Check the quality of these files using FASTQC.
```
1. Check the quality score of sample1. Interpret the report.

2. Check the quality score of sample2. Interpret the report.
```


## Read trimming and filtering (Trimmomatic)
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

Things to do in Trimming:
* Quality trimming
* Adapter trimming

### 🌟 Hands-on: Trimming raw data

__Trimmomatic for paired-end (PE) reads__

The flow of reads when trimming with Trimmomatic Paired End mode

![](https://img-blog.csdnimg.cn/img_convert/cc0cdfae795b9c58631f7d954e34441f.png)


```bash
# Trimming reads with Trimmomatic. Notes that the adapter's at ILLUMINACLIP must be specified (try to figure out this by look into the tool's adapters dir)
trimmomatic PE \
-phred33 \
-threads 56 \
-trimlog $p_trim/S11.log \
-summary $p_trim/S11_sum.log \
$p_raw/S11_L001_R1_001.fastq.gz \
$p_raw/S11_L001_R2_001.fastq.gz \
$p_trim/S11_L001_R1_trimmed_paired.fastq.gz \
$p_trim/S11_L001_R1_trimmed_unpaired.fastq.gz \
$p_trim/S11_L001_R2_trimmed_paired.fastq.gz \
$p_trim/S11_L001_R2_trimmed_unpaired.fastq.gz \
ILLUMINACLIP:/mnt/rdisk/duydao/PROJECT/DNASEQ_BRCA12/ref/adapters/TruSeq3-PE-2.fa:2:30:10:8:3:true \
HEADCROP:5 \
CROP:140 \
LEADING:3 \
TRAILING:10 \
SLIDINGWINDOW:4:28 \
MINLEN:36
```
Remove unpaired reads
```bash
cd $p_trim

rm *_unpaired*
```

Check QC again after Trim:
```bash
fastqc S11_L001_R1_trimmed_paired.fastq.gz -o qc_checked/
fastqc S11_L001_R2_trimmed_paired.fastq.gz -o qc_checked/
```


# 3. Alignment: Mapping reads to reference genome

First, Index the reference genome (We have done before)

```bash
bwa index -a bwtsw hg38.fa
```


## Introduce to SAM/BAM file
### Align with BWA mem

```bash
# Align with BWA mem, using 4 threads
bwa mem -t 4 \
-R '@RG\tID:rg1\tSM:S11\tPL:illumina\tLB:lib1\tPU:M07220:1:TGCAGCTA+CCTAGAGT' \
$p_ref/hg38.fa \
$p_trim/S11_L001_R1_trimmed_paired.fastq.gz \
$p_trim/S11_L001_R2_trimmed_paired.fastq.gz > $p_align/S11_aln.sam
```
### Convert SAM to BAM
BAM is a binary compressed of SAM file, which is less heavier than SAM.
```bash
samtools view -Sb S11_aln.sam > S11_aln.bam

# Now we have the compress bam file, we can remove unused sam file.
rm S11_aln.sam
```

__Samtools bitflags__

The "bitflag" field is a 16-bit integer that encodes various properties of a mapped read.

Use this link to explore more about this:
https://broadinstitute.github.io/picard/explain-flags.html

```bash
# not primary alignment (0x100) & supplementary alignment (0x800)
samtools view -F 0x900 S11_aln.bam

# not primary alignment (0x100)
samtools view -f 0x100 S11_aln.bam
```

__Flagstat__

The flagstat function of SAMtools provides a summary of the number of records corresponding to each of the bit flags.

```bash
samtools flagstat S11_aln.bam

#
1629692 + 0 in total (QC-passed reads + QC-failed reads)
1629194 + 0 primary
0 + 0 secondary
498 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1618417 + 0 mapped (99.31% : N/A)
1617919 + 0 primary mapped (99.31% : N/A)
1629194 + 0 paired in sequencing
814597 + 0 read1
814597 + 0 read2
1603792 + 0 properly paired (98.44% : N/A)
1617768 + 0 with itself and mate mapped
151 + 0 singletons (0.01% : N/A)
108 + 0 with mate mapped to a different chr
53 + 0 with mate mapped to a different chr (mapQ>=5)
```

---
# 4. Mapped reads post-processing
## 🌟 [Hands-on]: Use GATK to analyze data

## Sorting and marking duplicate, inxdexing the BAM file

### SortSam
```bash
gatk SortSam \
--INPUT S11_aln.bam \
--OUTPUT S11_aln_sorted.bam \
--SORT_ORDER coordinate
```
    
### MarkDuplicates

Detect duplicate using GATK.
```bash
gatk MarkDuplicates \
--INPUT S11_aln_sorted.bam \
--OUTPUT S11_aln_dedup.bam \
--METRICS_FILE S11_MarkDup.metrics
```

```bash
# Full BAM file (with duplicate reads are marked)
samtools view S11_aln_dedup.bam

# BAM file with duplicate reads only
samtools view -f 0x400 S11_aln_dedup.bam
samtools view -c -f 0x400 S11_aln_dedup.bam
```

__PCR duplication & Optical Duplication__
```bash
# Marked the read as an optical duplicated (DT:Z:SQ) & PCR duplicated (DT:Z:LB)
gatk MarkDuplicates \
--INPUT S11_aln_sorted.bam \
--OUTPUT S11_aln_dedup.bam \
--METRICS_FILE S11_MarkDup_PCR_OPT.metrics \
--TAGGING_POLICY All

# View Optical duplicates
samtools view -f 0x400 S11_aln_dedup.bam | \
grep DT:Z:SQ | less -S

# View PCR duplicates
samtools view -f 0x400 S11_aln_dedup.bam | \
grep DT:Z:LB | less -S
```

### Collect Insert Size
Using GATK modules
```bash
# Insert size
gatk CollectInsertSizeMetrics \
-I S11_aln.bam \
-O S11_aln_InsertSize_metrics.txt \
-H S11_aln_InsertSize_hist.pdf
```

# 5. Alignment Data: Quality Control
## Coverage Analysis
### Case 1: Keep duplicated reads
To check the coverage of BRCA1 and BRCA2 regions, we must have their location which stored in .bed file.

Download the .bed file that contains BRCA1/2 locations.
https://www.neb.com/en/tools-and-resources/usage-guidelines/grch38hg38-bed-files-for-the-nebnext-direct-brca1-brca2-panel

```bash
# Calculate coverage using BED tools
bedtools coverage \
-hist \
-a $p_ref/bedfile/brca12.bed \
-b $p_align/S11_aln_sorted.bam > S11.bed.cov

# Extract all coverage
grep ^all S11.bed.cov > S11.all.cov
```

### Case 2: Remove duplicated reads
Remove duplicate
```bash
gatk MarkDuplicates \
--INPUT S11_aln_sorted.bam \
--OUTPUT S11_aln_remove_dup.bam \
--METRICS_FILE S11_remove_dup.metrics \
--REMOVE_DUPLICATES true
```

Calculate coverage using BED tools
```bash
bedtools coverage \
-hist \
-a $p_ref/bedfile/brca12.bed \
-b $p_align/S11_aln_remove_dup.bam > S11_rmdup_bed.cov

# Extract all coverage
grep ^all NIST.bed.cov > NIST.all.cov
grep ^all S11_rmdup_bed.cov > S11_rmdup.all.cov
```

## A coverage plot in R
Finally, plot the coverage.

### Case 1
```R
cover <- read.table("S11.all.cov")
cov_cumul <- 1-cumsum(cover[,5])
plot(cover[1:200, 2], cov_cumul[1:200], type='l',
xlab="Depth",
ylab="Fraction of capture target bases >= depth",
ylim=c(0,1.0),
col="red",
main="Target Region Coverage")
abline(v = 20, col = "gray60")
abline(v = 100, col = "gray60")
abline(v = 300, col = "gray60")
abline(v = 500, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,100,300,500), labels=c(20,100,300,500))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))
dev.off()
```

### Case 2
```R
cover <- read.table("S11_rmdup.all.cov")
cov_cumul <- 1-cumsum(cover[,5])
plot(cover[1:200, 2], cov_cumul[1:200], type='l',
xlab="Depth",
ylab="Fraction of capture target bases >= depth",
ylim=c(0,1.0),
col="red",
main="Target Region Coverage")
abline(v = 20, col = "gray60")
abline(v = 100, col = "gray60")
abline(v = 300, col = "gray60")
abline(v = 500, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,100,300,500), labels=c(20,100,300,500))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))
dev.off()
```
---------------------------------

