# Practice with variant calling, filtering and annotation

> One can perform all steps without having to install anything through this notebook: [Variant_analysis in cloud](https://colab.research.google.com/drive/1m4pCPuewYHemnWwn31kOEDHt0rB5hOgk?usp=drive_link).

## Table of Contents

- [Practice with variant calling, filtering and annotation](#practice-with-variant-calling-filtering-and-annotation)
  - [Table of Contents](#table-of-contents)
  - [1. Reference genome preparation](#1-reference-genome-preparation)
  - [2. Variant calling with Haplotypecaller](#2-variant-calling-with-haplotypecaller)
  - [3. Select variants](#3-select-variants)
    - [3.1 SNPs](#31-snps)
    - [3.2 Short indels](#32-short-indels)
  - [4. Variant hard filtering](#4-variant-hard-filtering)
    - [4.1 SNPs](#41-snps)
    - [4.2 Indels](#42-indels)
  - [5. Variant annotation with Funcotator](#5-variant-annotation-with-funcotator)
    - [5.1 Germline pre-packaged hg38 data sources](#51-germline-pre-packaged-hg38-data-sources)
    - [5.2 SNPs annotation](#52-snps-annotation)
    - [5.3 Indels annotation](#53-indels-annotation)

## 1. Reference genome preparation

```bash
samtools faidx "chr13_chr17.fasta.gz"
gatk CreateSequenceDictionary -R "chr13_chr17.fasta.gz"
```

## 2. Variant calling with Haplotypecaller

```bash
gatk --java-options "-Xmx4G" HaplotypeCaller \
    --reference "ref/chr13_chr17.fasta.gz" \
    --input "bam/S11.unmrk.recal.bam" \
    --output "variant_calling/S11.haplotypecaller.vcf.gz" \
    --intervals "ref/BRCA_intervals.bed" \
    -stand-call-conf 30
```

## 3. Select variants

### 3.1 SNPs

```bash
gatk SelectVariants \
        -R "ref/chr13_chr17.fasta.gz" \
        -V "variant_calling/S11.haplotypecaller.vcf.gz" \
        --select-type SNP \
        -O "variant_calling/S11.haplotypecaller.snps.vcf.gz"
```

### 3.2 Short indels

```bash
gatk SelectVariants \
        -R "ref/chr13_chr17.fasta.gz" \
        -V "variant_calling/S11.haplotypecaller.vcf.gz" \
        --select-type INDEL \
        -O "variant_calling/S11.haplotypecaller.indels.vcf.gz"
```

## 4. Variant hard filtering

### 4.1 SNPs

```bash
gatk VariantFiltration \
        -R "ref/chr13_chr17.fasta.gz" \
        -V "variant_calling/S11.haplotypecaller.snps.vcf.gz" \
        -O "variant_calling/S11.haplotypecaller.snps.filtered.vcf.gz" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 3.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
```

### 4.2 Indels

```bash
gatk VariantFiltration \
        -R "ref/chr13_chr17.fasta.gz" \
        -V "variant_calling/S11.haplotypecaller.indels.vcf.gz" \
        -O "variant_calling/S11.haplotypecaller.indels.filtered.vcf.gz" \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0"
```

## 5. Variant annotation with Funcotator

### 5.1 Germline pre-packaged hg38 data sources

```bash
gatk FuncotatorDataSourceDownloader --germline --validate-integrity --hg38 --extract-after-download
```

### 5.2 SNPs annotation

```bash
gatk --java-options "-Xmx4G" Funcotator \
        -V "variant_calling/S11.haplotypecaller.snps.filtered.vcf.gz" \
        -R "ref/chr13_chr17.fasta.gz" \
        -O "annotation/S11.haplotypecaller.snps.filtered.ann.vcf.gz" \
        --output-file-format VCF \
        --data-sources-path "path/to/funcotator_dataSources.v1.8.hg38.20230908g" \
        --ref-version hg38

# Parsing result into tab-delimited file
gatk VariantsToTable \
    -V "annotation/S11.haplotypecaller.snps.filtered.ann.vcf.gz" -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O "annotation/S11.haplotypecaller.snps.filtered.ann.tsv"
```

### 5.3 Indels annotation

```bash
gatk --java-options "-Xmx4G" Funcotator \
        -V "variant_calling/S11.haplotypecaller.indels.filtered.vcf.gz" \
        -R "ref/chr13_chr17.fasta.gz" \
        -O "annotation/S11.haplotypecaller.indels.filtered.ann.vcf.gz" \
        --output-file-format VCF \
        --data-sources-path "path/to/funcotator_dataSources.v1.8.hg38.20230908g" \
        --ref-version hg38

# Parsing result into tab-delimited file
gatk VariantsToTable \
    -V "annotation/S11.haplotypecaller.snps.indels.ann.vcf.gz" -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O "annotation/S11.haplotypecaller.indels.filtered.ann.tsv"
```
