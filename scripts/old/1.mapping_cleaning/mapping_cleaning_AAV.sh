#!/bin/bash

# Version 1.0
# This script was last modified on 13/10/2021
# Script Goal: AAV mapping cleaning
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After mapping with bwa-mem, filtering the output
# Remove reads due to technical and biological biases
# Identify regions too similar to mouse genome
# Identify regions containing repeat elements present in mouse genome
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file of interest in quotes
# 3. Write the name of the vector genome reference fasta file of interest in quotes
# 4. Write the name of the mouse genome reference fasta file of interest in quotes
#
# Output: ???
# 
# Concept:
# The default method followed is normalization via division by the sum of sequences in a given sample
# 
# Note: ???


#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/mapping_cleaning/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/1.mapping_cleaning/output

# Give the file name of the original vector genome
vector_genome="${workdir}/0.preparing_input_files/vector_genome.fasta" #------------------------> !CHANGE IT!#
mouse_genome="${workdir}/0.preparing_input_files/mouse_genome.fa" #------------------------> !CHANGE IT!#

# Give the file name of the AAV mapping BAM file
mapping_file="${workdir}/0.preparing_input_files/149_mapping_AAV.bam"

# Give the file name of repeat masker result file
repeat_elements="${workdir}/0.preparing_input_files/AAV_RepeatMasker.fa.out"

# Set the number of threads to use
nb_threads=100

# Set the name of the sample
sample=154 #------------------------> !CHANGE IT!#

# Create a folder to store all the results for this sample
mkdir -p $workdir/1.mapping_cleaning/output/vector/${sample}



#----------------- END: No changes are needed below this line --------------------#



##################################################################################
######                             Main Script                              ######
##################################################################################

#----------------- Filter out reads due to technical biases (unmapped, duplicates, multi-mapped) --------------------#

# Filter out alignments which correspond to the FLAG option
## Tool: samtools view
## Options:
### -b: output in the BAM format
### -F: filter out alignments with all bits in INT present in the FLAG field
### -@: number of threads

### FLAG: 4 -> remove unmapped reads from the original AAV mapping BAM file
../bin/samtools view -@ $nb_threads -b -F 4 $mapping_file > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.bam

# Identify and remove duplicate reads (PCR and optical duplicates) from the AAV mapping BAM file
## Tool: sambamba markdup
## Options:
### -r: remove duplicates instead of just marking them
### -t: number of threads
### -p: show progressbar in the terminal
../bin/sambamba markdup -r -t $nb_threads -p $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.bam $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.bam

### FLAG: 256 -> remove secondary alignments from the AAV mapping BAM file
../bin/samtools view -@ $nb_threads -b -F 256 $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.bam > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.bam 

### FLAG: 2048 -> remove supplementary alignments from the AAV mapping BAM file
../bin/samtools view -@ $nb_threads -b -F 2048 $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.bam > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.bam

### Options: -q 1 -> remove reads with a mapping quality equal to zero (lower than 1)
../bin/samtools view -b -q 1 $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.bam > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam



#----------------- Filter out reads which are associated with regions too similar to the mouse genome --------------------#

# Align the AAV genome reference against the mouse reference to identify is there is similar regions between the two genomes
## Tool: blastn
## Options:
### -query: AAV vector genome reference
### -subject: mouse vector genome reference
### -out: name of the output file
### -outfmt: output format as a table with headers (and not in the form of alignments)
../bin/blastn -query $vector_genome -subject $mouse_genome -out $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.txt -outfmt 6

# Adding an header of the output file
sed -i -e '1i\qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.txt

# Prepare the BED file with only three columns: chr, start and end of each similar regions
sed '1d' $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.txt | awk '{print $1 "\t" $7 "\t" $8}' > $workdir/1.mapping_cleaning/output/vector_similar_regions_to_remove.bed

# Remove reads where overlap is found
## Tool: bedtools intersect
## Options:
### -v: keep only reads in the AAV mapping file that have no overlap in the regions identified in the BED file
### -abam: each read (BAM alignment) in the AAV mapping file is compared to the regions in the BED file in search of overlaps
### -b: to specify the following BED file
### MODIFICATION! blastn-short
#../bin/bedtools intersect -v -abam $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam -b $workdir/1.mapping_cleaning/output/vector_similar_regions_to_remove.bed > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.bam
../bin/bedtools intersect -v -abam $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam -b $workdir/1.mapping_cleaning/output/vector_similar_regions_to_remove_short.bed > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam



#----------------- Filter out reads which are associated with repeat element identified in the mouse genome --------------------#

# Take the result from Repeat Masker analysis on the annotation of the repeats that are present in the query sequence
## Tool: repeat masker
## Species: mouse genome

# Prepare the BED file with only three columns: chr, start and end of each repeat regions that correspond to the repeat elements in the mouse genome
sed '1,3d' $repeat_elements | awk '{print $6 "\t" $7}' | awk 'BEGIN{FS=OFS="\t"} {print "chr1", $0}' > $workdir/1.mapping_cleaning/output/vector_repeat_regions_to_remove.bed

# Remove reads where overlap is found
## Tool: bedtools intersect
## Options:
### -v: keep only reads in the AAV mapping file that have no overlap in the regions identified in the BED file
### -abam: each read (BAM alignment) in the AAV mapping file is compared to the regions in the BED file in search of overlaps
### -b: to specify the following BED file
../bin/bedtools intersect -v -abam $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.bam -b $workdir/1.mapping_cleaning/output/vector_repeat_regions_to_remove.bed > $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam 


#------------------------------------------- END --------------------------------------------------#