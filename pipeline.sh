#!/bin/bash

# AUTHOR: Ilja Lomkov
# STUDENT NUMBER: 1810106
# GITHUB LINK: https://github.com/Ilja-Lom
# DATE: 08/04/2022
# DESCRIPTION: A BASH script pipeline for end-to-end processing of FASTQ reads

# NOTES **************************************************
## This BASH script must be located within the HOME directory
## This script assumes that all the prerequisites are installed

# PREREQUISITES ******************************************
## ---- Anaconda ----
### samtools
### bwa
### freebayes
### picard
### bedtools
### trimmomatic
### fastqc
### vcflib

## ---- ANNOVAR ----
### refGene
### snp129

## 30GB of storage

# SETTING UP THE DIRECTORIES *****************************
## ensure current directory is home (cd ~/)
echo "GENERATING REQUIRED DIRECTORIES"
# creating a directory to store the assignment workings
mkdir assignment

# directories to store data, results, and logs
mkdir assignment/{data,results,logs}

# sub-directories within the 'data' directory
mkdir assignment/data/{reference,bed,fastq_untrimmed,fastq_trimmed,aligned_data,variants}

# sub-directories within the 'results' directory
mkdir assignment/results/{fastqc_untrimmed,fastqc_trimmed,post_alignment_metrics,annovar}

# sub-directories within the 'logs' directory
mkdir assignment/logs/{fastqc_untrimmed_summaries,fastqc_trimmed_summaries,annovar_output}

# DOWNLOADING DATA ***************************************

# downloading the FASTQ READ 1
echo "DOWNLOADING READ 1"
wget -O ~/assignment/data/fastq_untrimmed/NGS0001_R1.fastq.gz https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
gunzip ~/assignment/data/fastq_untrimmed/NGS0001_R1.fastq.gz

# downloading the FASTQ READ 2
echo "DOWNLOADING READ 2"
wget -O ~/assignment/data/fastq_untrimmed/NGS0001_R2.fastq.gz https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
gunzip ~/assignment/data/fastq_untrimmed/NGS0001_R2.fastq.gz

# downloading the reference file
echo "DOWNLOADING REFERENCE"
wget -O ~/assignment/data/reference/hg19.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# downloading the annotation BED
echo "DOWNLOADING BED FILE"
wget -O ~/assignment/data/bed/annotation.bed https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# ASSESSING THE FASTQC QUALITY OF UNTRIMMED FASTQ READS *******************

# running FASTQC on all the UNTRIMMED FASTQ files
echo "GENERATING FASTQC OUTPUT ON UNTRIMMED > ~/assignment/logs/fastqc_untrimmed_summaries"
fastqc -t 4 ~/assignment/data/fastq_untrimmed/*fastq -o ~/assignment/logs/fastqc_untrimmed_summaries
echo "FASTQC REPORT SUCCESSFULLY GENERATED"
# TRIMMING THE FASTQ READS ************************************************
## using Trimmomatic

# ---- INPUTS ----
## ~/assignment/data/fastq_untrimmed/NGS0001.R1.fastq.qz
## ~/assignment/data/fastq_untrimmed/NGS0001.R2.fastq.qz

# ---- OUTPUT ----
## ~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads_1P (.fastq)
## ~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads_2P (.fastq)

echo "TRIMMOMATIC TRIMMING FASTQ READS"
trimmomatic PE \
-threads 4 \
-phred33 \
~/assignment/data/fastq_untrimmed/NGS0001_R1.fastq \
~/assignment/data/fastq_untrimmed/NGS0001_R2.fastq \
-baseout ~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads \
ILLUMINACLIP:/home/ubuntu/miniconda/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10:2:True \
LEADING:3 TRAILING:3 MINLEN:36
echo "TRIMMING SUCCESSFUL"
# ASSESSING THE FASTQC QUALITY OF TRIMMED FASTQ READS *********************

# running FASTQC on all the TRIMMED FASTQ files
echo "GENERATING FASTQC OUTPUT ON TRIMMED > ~/assignment/logs/fastqc_trimmed_summaries"
fastqc -t 4 ~/assignment/data/fastq_trimmed/*P -o ~/assignment/logs/fastqc_trimmed_summaries
echo "FASTQC REPORT SUCCESSFULLY GENERATED"

# ZIPPING THE UNTRIMMED FASTQ FILES TO SAVE SPACE

# ---- INPUT ----
## ~/assignment/data/fastq_untrimmed/NGS0001_R1.fastq
## ~/assignment/data/fastq_untrimmed/NGS0001_R2.fastq

# ---- OUTPUT ----
## ~/assignment/data/fastq_untrimmed/NGS0001_R1.fastq.gz
## ~/assignment/data/fastq_untrimmed/NGS0001_R2.fastq.gz

echo "COMPRESSING UNTRIMMED FASTQ"
gzip ~/assignment/data/fastq_untrimmed/NGS0001_R1.fastq
gzip ~/assignment/data/fastq_untrimmed/NGS0001_R2.fastq
echo "COMPRESSION SUCCESSFUL"

# INDEXING THE REFERENCE SEQUENCE USING BWA INDEX *************************

# ---- INPUT ----
## ~/assignment/data/reference/hg19.fa.gz

# ---- OUTPUT ----
## NULL

echo "BWA INDEXING THE REFERENCE"
bwa index ~/assignment/data/reference/hg19.fa.gz
echo "INDEXING SUCCESSFUL"

# ALIGNMENT ***************************************************************
## aligning the trimmed reads to the indexed reference sequence

# ---- INPUT ----
## ~/assignment/data/reference/hg19.fa.gz
## ~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads_1P.fastq
## ~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads_2P.fastq

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned.sam

echo "BWA ALIGNING TRIMMED READS TO REFERENCE"
bwa mem -t 4 -v 1 -R \
'@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 \
~/assignment/data/reference/hg19.fa.gz \
~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads_1P \
~/assignment/data/fastq_trimmed/NGS0001_trimmed_reads_2P \
> ~/assignment/data/aligned_data/NGS0001_aligned_sam.sam
echo "ALIGNMENT SUCCESSFUL"

# CONVERTING SAM TO BAM ***************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned.sam

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned.bam

echo "PERFORMING SAMTOOLS SAM>BAM: ~/assignment/data/aligned_data/NGS0001_aligned_sam.sam > ~/assignment/data/aligned_data/NGS0001_aligned.bam"
samtools view -h -b ~/assignment/data/aligned_data/NGS0001_aligned_sam.sam > ~/assignment/data/aligned_data/NGS0001_aligned.bam
echo "SAMTOOLS SAM>BAM SUCCESSFUL"
# removing the SAM file to save space
echo "REMOVING REDUNDANT: ~/assignment/data/aligned_data/NGS0001_aligned_sam.sam"
rm ~/assignment/data/aligned_data/NGS0001_aligned_sam.sam
echo "REMOVAL SUCCESSFUL"

# SORTING BAM FILE ********************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned.bam

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam

echo "PERFORMING SAMTOOLS SORT: ~/assignment/data/aligned_data/NGS0001_aligned.bam > ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam"
samtools sort ~/assignment/data/aligned_data/NGS0001_aligned.bam > ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam
echo "SAMTOOLS SORT COMPLETED"

echo "REMOVING REDUNDANT: ~/assignment/data/aligned_data/NGS0001_aligned.bam"
rm ~/assignment/data/aligned_data/NGS0001_aligned.bam
echo "REMOVAL SUCCESSFUL"

# INDEXING SORTED BAM FILE ************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam.bai

echo "RUNNING SAMTOOLS INDEX: ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam > ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam.bai"
samtools index ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam > ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam.bai
echo "SAMTOOLS INDEX COMPLETED SUCCESSFULLY"

# POST-ALIGNMENT QUALITY CONTROL ******************************************
# MARKING DUPLICATES ******************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam 

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam
## ~/assignment/results/post_alignment_metrics/marked_duplicate_metrics.txt

# marking the duplicates in the aligned and sorted reads
echo "RUNNING PICARD MARKDUPLICATES"
picard MarkDuplicates I=/home/ubuntu/assignment/data/aligned_data/NGS0001_aligned_sorted.bam O=/home/ubuntu/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam M=/home/ubuntu/assignment/results/post_alignment_metrics/marked_duplicate_metrics.txt
echo "PICARD MARKDUPLICATES COMPLETED SUCCESSFULLY"

# REMOVING DUPLICATES *****************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam.bai

# removing the duplicates marked by Picard
echo "REMOVING DUPLICATES"
samtools index ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam > ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam.bai
echo "DUPLICATES REMOVED SUCCESSFULLY"

# FILTERING BAM BY MAPPING QUALITY ****************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam

echo "SAMTOOLS VIEW: FILTERING BAM BY MAP QUALITY"
## NOTE: format is output first, then input
samtools view -F 1796 -q 20 -o ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam
echo "FILTERING COMPLETED SUCCESSFULLY"

# INDEXING FILTERED BAM ***************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam

# ---- OUTPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam.bai

echo "SAMTOOLS INDEX: INDEXING FILTERED BAM > ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam"
samtools index ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam > ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam.bai
echo "INDEXING SUCCESSFUL > ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam.bai"

# INSERT ADDITIONAL POST-ALIGNMENT-METRIC TESTS BELOW ********************************

# FLAGSTATS

# samtools flagstat ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam -O ~/assignment/results/post_alignment_metrics/flagstats.tsv

# IDXSTATS

# samtools idxstats ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam > ~/assignment/results/post_alignment_metrics/idxstats.tsv

# DEPTH OF COVERAGE *************************************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam
## ~/assignment/data/reference/hg19.fa

# ---- OUTPUT ----
## ~/assignment/results/post_alignment_metrics/depth_of_coverage.png

echo "RUNNING BEDTOOLS GENOMECOV"
bedtools genomecov -ibam ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam -g ~/assignment/data/reference/hg19.fa > ~/assignment/results/post_alignment_metrics/depth_of_coverage.txt
echo "COMPLETED SUCCESSFULLY > depth_of_coverage.png"

# INSERT SIZE *****************************************************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam

# ---- OUTPUT ----
## ~/assignment/results/post_alignment_metrics/insert_size_metrics.txt
## ~/assignment/results/post_alignment_metrics/insert_size_histogram.pdf

echo "RUNNING PICARD COLLECT_INSERT_SIZE_METRICS"
picard CollectInsertSizeMetrics \
I=~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam \
O=~/assignment/results/post_alignment_metrics/insert_size_metrics.txt \
H=~/assignment/results/post_alignment_metrics/insert_size_histogram.pdf \
M=0.5
echo "COMPLETED SUCCESSFULLY > ~/assignment/results/post_alignment_metrics/insert_size_metrics.txt + ~/assignment/results/post_alignment_metrics/insert_size_histogram.pdf"

# REMOVING UNNECESSARY FILES ************************************************
## saving space

echo "REMOVING UNNECESSARY FILES"
echo "REMOVING > ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam"
rm ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam
echo "REMOVING > ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam.bai"
rm ~/assignment/data/aligned_data/NGS0001_aligned_sorted.bam.bai
echo "REMOVING > ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam"
rm ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam
echo "REMOVING > ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam.bai"
rm ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked.bam.bai

## NOTE: this saves 3GB of storage

# VARIANT CALLING AND FILTERING *******************************************
# *************************************************************************

# UNZIP REFERENCE *********************************************************

# ---- INPUT ----
## ~/assignment/data/reference/hg19.fa.gz

# ---- OUTPUT ----
## ~/assignment/data/reference/hg19.fa

echo "UNZIPPING REFERENCE > ~/assignment/data/reference/hg19.fa.gz"
zcat ~/assignment/data/reference/hg19.fa.gz > ~/assignment/data/reference/hg19.fa
echo "UNZIPPING SUCCESSFUL > ~/assignment/data/reference/hg19.fa"

# INDEXING THE UNZIPPED REFERENCE *****************************************

# ---- INPUT ----
## ~/assignment/data/reference/hg19.fa

# ---- OUTPUT ----
## ~/assignment/data/reference/hg19.fa.fai

echo "SAMTOOLS FAIDX: INDEXING REFERENCE > ~/assignment/data/reference/hg19.fa"
samtools faidx ~/assignment/data/reference/hg19.fa > ~/assignment/data/reference/hg19.fa.fai
echo "INDEXING SUCCESSFUL > ~/assignment/data/reference/hg19.fa.fai"

# FREEBAYES VARIANT CALLING ***********************************************

# ---- INPUT ----
## ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam
## ~/assignment/data/reference/hg19.fa

# ---- OUTPUT ----
## ~/assignment/data/variants/NGS0001_vcf.vcf

echo "FREEBAYES: VARIANT CALLING"
freebayes --bam ~/assignment/data/aligned_data/NGS0001_aligned_sorted_marked_filtered.bam --fasta-reference ~/assignment/data/reference/hg19.fa --vcf ~/assignment/data/variants/NGS0001_vcf.vcf
echo "FREEBAYES VARIANT CALLING SUCCESSFUL > ~/assignment/data/variants/NGS0001_vcf.vcf"

# COMPRESSING FREEBAYES VCF OUTPUT ****************************************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_vcf.vcf

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_vcf.vcf.gz

echo "COMPRESSING FREEBAYES VCF OUTPUT > ~/assignment/data/variants/NGS0001_vcf.vcf"
bgzip ~/assignment/data/variants/NGS0001_vcf.vcf
echo "COMPRESSION SUCCESSFUL > ~/assignment/data/variants/NGS0001_vcf.vcf.gz"

# NOTE: the input file is removed after compressing, no need to input an output file.
## SOURCE: http://www.htslib.org/doc/bgzip.html

# INDEXING THE VARIANT VCF FILE *******************************************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_vcf.vcf.gz

# ---- OUTPUT ----
## ~/assignment/data/variants/NGS0001_vcf.vcf.gz.tbi

echo "TABIX: INDEXING VARIANT VCF.GZ FILE > ~/assignment/data/variants/NGS0001_vcf.vcf.gz"
tabix -p vcf ~/assignment/data/variants/NGS0001_vcf.vcf.gz
echo "TABIX: INDEXING OF VARIANT VCF.GZ SUCCESSFUL > ~/assignment/data/variants/NGS0001_vcf.vcf.gz.tbi"

# FILTERING THE VCF *******************************************************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_vcf.vcf.gz

# ---- OUTPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf

echo "VCFFILTER: FILTERING VARIANT VCF.GZ > ~/assignment/data/variants/NGS0001_vcf.vcf.gz"
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/assignment/data/variants/NGS0001_vcf.vcf.gz > ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf
echo "VCFFILTER: VARIANT FILTERING SUCCESSFUL > ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf"

# FILTERING BY BED FILE ****************************************************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf
## ~/assignment/data/bed/annotation.bed

# ---- OUTPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf

echo "BEDTOOLS: FILTERING VCF BY BED FILE > ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf + ~/assignment/data/bed/annotation.bed"
bedtools intersect -header -wa -a ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf -b ~/assignment/data/bed/annotation.bed > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf
echo "BEDTOOLS: FILTERING VCF BY BED SUCCESSFUL > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf"

# ZIPPING FILTERED **********************************************************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf

# ---- OUTPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz

echo "COMPRESSING BED FILTERED VARIANTS > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf"
bgzip ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf
echo "COMPRESSION SUCCESSFUL > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz"

## NOTE: ORIGINAL FILE REMOVED

# INDEXING THE FILTERED BED VCF *********************************************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz

# ---- OUTPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz.tbi

echo "TABIX: INDEXING FILTERED BED VARIANTS > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz"
tabix -p vcf ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz
echo "TABIX: INDEXING SUCCESSFUL > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz.tbi"

# REMOVING UNNECESARY FILES - SAVING STORAGGE SPACE *************************
echo "REMOVING UNNECESSARY FILES"

echo "REMOVING > ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf"
rm ~/assignment/data/variants/NGS0001_filtered_final_vcf.vcf
echo "REMOVING > ~/assignment/data/variants/NGS0001_vcf.vcf.gz"
rm ~/assignment/data/variants/NGS0001_vcf.vcf.gz
echo "REMOVING > ~/assignment/data/variants/NGS0001_filtered_vcf.vcf.gz.tbi"
rm ~/assignment/data/variants/NGS0001_vcf.vcf.gz.tbi
echo "REMOVAL SUCCESSFUL"

## NOTE: this frees 19GB of data

# ANNOVAR ANNOTATION ********************************************************
# ***************************************************************************

echo "ANNOVAR ANNOTATION"


# CONVERTING COMPRESSED VCF FILE INTO ACCEPTED FORMAT ***********************

# ---- INPUT ----
## ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz

# ---- OUTPUT ----
## ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.avinput

echo "ANNOVAR: CONVERTING VCF TO ANNOVAR > ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz"
~/annovar/convert2annovar.pl -format vcf4 ~/assignment/data/variants/NGS0001_filtered_final_bed_vcf.vcf.gz > ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.avinput
echo "ANNOVAR: CONVERSION SUCCESSFUL > ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.avinput"

# RUNNING ANNOVAR ***********************************************************

# ---- INPUT ----
## ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.avinput

# ---- OUTPUT ----
## ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.hg19_multianno.csv

echo "ANNOVAR: RUNNING ANNOVAR"
~/annovar/table_annovar.pl ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.avinput annovar/humandb/ -buildver hg19 -out ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf -remove -protocol refGene,snp129 -operation g,f -nastring . -csvout
echo "ANNOVAR: ANNOVAR PROCESSING SUCCESSFUL > ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.hg19_multianno.csv"

# FILTERING ANNOTATED VARIANTS BY EXONIC AND UNSEEN IN dbSNP129 *******************

## Due to storage space limitations, an outdated version of dbSNP is used (dbSNP129)

# ---- INPUT ----
## ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.hg19_multianno.csv

# ---- OUTPUT ----
## ~/assignment/results/annovar/NGS0001_filtered_annovar.csv

# use awk to filter for exonic variants in the 6th column
# indicate csv, comma delimited, by stating -F ","
# input csv has quotes around exonic ("exonic"), a backslash indicates it is part of the string
# print command will return ALL columns

echo "FILTERING ANNOVAR ANNOTATION: EXONIC, UNSEEN IN DBSNP"
awk -F "," '{if( ($6 == "\"exonic\"") && ($11 == ".") ) print}' ~/assignment/results/annovar/NGS0001_filtered_final_bed_vcf.hg19_multianno.csv > ~/assignment/results/annovar/NGS0001_filtered_annovar.csv
echo "FILTERING ANNOVAR SUCCESSFUL > ~/assignment/results/annovar/NGS0001_filtered_annovar.csv"

echo "PIPELINE COMPLETED"
# END **************************************************************************




































