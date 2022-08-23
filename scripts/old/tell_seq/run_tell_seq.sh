#!/bin/bash

# This script allow to launch Tell-Read (mapping to both genomes) pipeline on the samples specified in a file: 'sample_list.txt'
# Attention! Take 2 To of free space when running one sample at a time


#------------------------------ Set up parameters ------------------------------------#

# Specify the working directory
workdir="/data/kevin/Tellseq"

# Specify the genome reference
genome_AAV="AAV_hGPE_HBB2_hG6PCwt_bGH"
genome_mouse="GRCm38"


#------------------------------ Prepare docker image (only one time!) ------------------------------------#

# # Load the docker image
# docker load -i docker-tellread
# docker load -i docker-tellsort

# # Check image docker is loaded
# docker images


#------------------------------ Prepare genome reference directory (only one time!) ------------------------------------#

## Create directories and put genome reference fasta inside
#mkdir genomes
#cd genomes/
#mkdir GRC38
#mkdir AAV_ref

## This script generates the genome's indexes and bed files for each genome reference (for AAV vector)
#$workdir/tellread-release/generateGenomeIndexBed.sh AAV_hGPE_HBB2_hG6PCwt_bGH.fasta

## This script generates the genome's indexes and bed files for each genome reference (for mouse)
#$workdir/tellread-release/generateGenomeIndexBed.sh GRCm38_68.fa



#--- START OF THE TELL-SEQ ANALYSIS ---#

# Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"


#------------------------------ Prepare files for Tell-Seq ------------------------------------#

# Extract 8 files from the compressed folder


# Identify R1, R2, I1 and I2 #

## Look at the file & rename it
zcat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz | head ## index 2 (8-base)
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I2.fastq.gz # first lane
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I2.fastq.gz # second lane

zcat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz | head ## index 1 (18-base)
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz # first lane
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I1.fastq.gz # second lane

zcat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R3.fastq.gz | head ## read 2
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R3.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz # first lane
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R3.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R2.fastq.gz # second lane

## Merge two fastq.gz files #
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I2.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I2.fastq.gz
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I1.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I1.fastq.gz
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R2.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R2.fastq.gz
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R1.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R1.fastq.gz

## Clean files with one lane
rm $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.*
rm $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.*
rm $workdir/GC1*.fastq.gz



#------------------------------ Run Tell-Read pipeline on fastq raw data ------------------------------#

# Subsampling & performance analysis on mouse genome
$workdir/tellread-release/run_tellread_fq.sh \
-i1 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I1.fastq.gz \
-i2 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I2.fastq.gz \
-r1 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R1.fastq.gz \
-r2 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R2.fastq.gz \
-o $workdir/${sample}/demultiplex_mouse \
-f $workdir/genomes/ \
-s T500 \
-g ${genome_mouse}

## CLEAN!
rm $workdir/${sample}/demultiplex_mouse/1_demult/Raw/demultiplex_mouse_*
# Remove files without zip extension
find $workdir/${sample}/demultiplex_mouse/1_demult/fastqc -type f ! -name '*.zip' -delete
# Compress large files
pigz $workdir/${sample}/demultiplex_mouse/2_barcode_indiv/demultiplex_mouse_T500_all_barcode_seq_count.txt
pigz $workdir/${sample}/demultiplex_mouse/2_barcode_indiv/demultiplex_mouse_T500_correct_barcode_seq_count.txt



# Subsampling & performance analysis on AAV reference sequence
$workdir/tellread-release/run_tellread_fq.sh \
-i1 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I1.fastq.gz \
-i2 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I2.fastq.gz \
-r1 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R1.fastq.gz \
-r2 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R2.fastq.gz \
-o $workdir/${sample}/demultiplex_AAV \
-f $workdir/genomes/ \
-s T500 \
-g ${genome_AAV}

## CLEAN!
rm $workdir/${sample}/demultiplex_AAV/1_demult/Raw/demultiplex_AAV_*
# Remove files without zip extension
find $workdir/${sample}/demultiplex_AAV/1_demult/fastqc -type f ! -name '*.zip' -delete
# Compress large files
pigz $workdir/${sample}/demultiplex_AAV/2_barcode_indiv/demultiplex_AAV_T500_all_barcode_seq_count.txt
pigz $workdir/${sample}/demultiplex_AAV/2_barcode_indiv/demultiplex_AAV_T500_correct_barcode_seq_count.txt



#------------------------------------------- Prepare genome VCF directory -----------------------------------------#

# Here, no reference VFC files are needed!



#---------------------------------- Run Tell-Sort pipeline on linked-read fastq data ----------------------------------#

# Align with the mouse genome
$workdir/tellsort-release/run_tellsort.sh \
-r1 $workdir/${sample}/demultiplex_mouse/Full/demultiplex_mouse_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
-r2 $workdir/${sample}/demultiplex_mouse/Full/demultiplex_mouse_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
-i1 $workdir/${sample}/demultiplex_mouse/Full/demultiplex_mouse_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
-r $workdir/genomes/${genome_mouse}/${genome_mouse}.fa \
-o $workdir/${sample}/align_mouse \
-p align_mouse \
-t 65

## CLEAN! Compress large files
pigz $workdir/${sample}/align_mouse/align_mouse_temp/align_mouse.r1.fastq
pigz $workdir/${sample}/align_mouse/align_mouse_temp/align_mouse.r2.fastq
pigz $workdir/${sample}/align_mouse/align_mouse_temp/molecule.tsv



# Align with the AAV reference vector sequence
$workdir/tellsort-release/run_tellsort.sh \
-r1 $workdir/${sample}/demultiplex_AAV/Full/demultiplex_AAV_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
-r2 $workdir/${sample}/demultiplex_AAV/Full/demultiplex_AAV_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
-i1 $workdir/${sample}/demultiplex_AAV/Full/demultiplex_AAV_I1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
-r $workdir/genomes/${genome_AAV}/${genome_AAV}.fasta \
-o $workdir/${sample}/align_AAV \
-p align_AAV \
-t 60

# CLEAN! Compress large files
pigz $workdir/${sample}/align_AAV/align_AAV_temp/molecule.filter.tsv 
pigz $workdir/${sample}/align_AAV/align_AAV_temp/molecule.tsv
pigz $workdir/${sample}/align_AAV/align_AAV_temp/align_AAV.r1.fastq
pigz $workdir/${sample}/align_AAV/align_AAV_temp/align_AAV.r2.fastq



#------------------------------- Clear up unnecessary files ------------------------------------#

#rm $workdir/${sample}/align_AAV/align_AAV_temp/align_AAV.r* 


done < sample_list.txt

#--- END OF THE TELL-SEQ ANALYSIS ---#