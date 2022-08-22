#!/bin/bash

# Version 1.0
# This script was last modified on 26/10/2021 (UPDATE! Remove similar regions with AAV genome)
# Script Goal: Mouse mapping cleaning
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After mapping with bwa-mem, filtering the output
# Remove reads due to technical and biological biases
# Identify regions too similar to AAV vector genome
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the mouse mapping BAM file of interest in quotes
# 3. Write the name of the vector genome reference fasta file of interest in quotes
# 4. Write the name of the mouse genome reference fasta file of interest in quotes
#
# Notes:
# ERROR! Problem with the BED file for similar regions between the two genomes: wrong order between start and end positions!!!
# 	Cause: match inverse regions (=> RESOLVED!)




#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/mapping_cleaning/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
sample=154 #------------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/1.mapping_cleaning/output

# Give the file name of the original vector genome
vector_genome="${workdir}/0.preparing_input_files/vector_genome.fasta" #------------------------> !CHANGE IT!#
mouse_genome="${workdir}/0.preparing_input_files/mouse_genome.fa" #------------------------> !CHANGE IT!# (replace with the new one!)

# Give the file name of the mouse mapping BAM file (ln -s if the mapping file is too big to be deplaced!)
# Path to file: /media/kevincheeseman/WLG/149/align_mouse/align_mouse_temp
# Path to send: /data/kevin/Tellseq/ISA/0.preparing_input_files/
# CMD: ln -s /media/kevincheeseman/WLG/149/align_mouse/align_mouse_temp/
mapping_file="${workdir}/0.preparing_input_files/${sample}_mapping_mouse.bam" #------------------------> !CHANGE IT!#

# Give the file name of repeat masker result file
repeat_elements="${workdir}/0.preparing_input_files/mm10_RepeatMasker.fa.out"

# Set the number of threads to use
nb_threads=10

# Create a folder to store all the results for this sample
mkdir -p $workdir/1.mapping_cleaning/output/mouse/${sample}


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
../bin/samtools view -@ $nb_threads -b -F 4 $mapping_file > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.bam

# Identify and remove duplicate reads (PCR and optical duplicates) from the AAV mapping BAM file
## Tool: sambamba markdup
## Options:
### -r: remove duplicates instead of just marking them
### -t: number of threads
### -p: show progressbar in the terminal
### --tmpdir: change temporary directory location while running sambamba markdup (UPDATE!)
../bin/sambamba markdup -r -t $nb_threads --overflow-list-size 6000000 --hash-table-size 1000000 --tmpdir $workdir/1.mapping_cleaning/output/mouse/${sample}/ -p $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.bam $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.bam

### FLAG: 256 -> remove secondary alignments from the AAV mapping BAM file
../bin/samtools view -@ $nb_threads -b -F 256 $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.bam > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.bam 

### FLAG: 2048 -> remove supplementary alignments from the AAV mapping BAM file
../bin/samtools view -@ $nb_threads -b -F 2048 $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.bam > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.bam

### Options: -q 1 -> remove reads with a mapping quality equal to zero (lower than 1)
../bin/samtools view -b -q 1 $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.bam > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam



#----------------- Filter out reads which are associated with regions too similar to the mouse genome --------------------#

# Align the AAV genome reference against the mouse reference to identify is there is similar regions between the two genomes
## Tool: blastn
## Options:
### -query: AAV vector genome reference
### -subject: mouse vector genome reference
### -out: name of the output file
### -outfmt: output format as a table with headers (and not in the form of alignments)
### ATTENTION! Problem with order
../bin/blastn -query $vector_genome -subject $mouse_genome -out $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.tmp -outfmt 6

# Tool: awk to solve the problem with order number
## Concept: if the 'start' position number is lower than the 'end' position number inverse the two numbers
awk -F "\t" '{if($2<$3) print $0, $4=$2, $5=$3, $6="OK"; else print $0, $4=$3, $5=$2, $6="ERROR"}' $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.tmp | awk '{print $1 "\t" $4 "\t" $5}' > $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.txt

# Adding an header of the output file
sed -i -e '1i\qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.txt

# Prepare the BED file with only three columns: chr, start and end of each similar regions
sed '1d' $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.txt | awk '{print $2 "\t" $9 "\t" $10}' > $workdir/1.mapping_cleaning/output/mouse_similar_regions_to_remove.bed

# Remove reads where overlap is found
## Tool: bedtools intersect
## Options:
### -v: keep only reads in the AAV mapping file that have no overlap in the regions identified in the BED file
### -abam: each read (BAM alignment) in the AAV mapping file is compared to the regions in the BED file in search of overlaps
### -b: to specify the following BED file
#../bin/bedtools intersect -v -abam $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam -b $workdir/1.mapping_cleaning/output/mouse_similar_regions_to_remove.bed > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.bam
../bin/bedtools intersect -v -abam $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam -b $workdir/1.mapping_cleaning/output/mouse_similar_regions_to_remove_short.bed > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam


# Clean
rm $workdir/1.mapping_cleaning/output/blast_output_vector_mouse.tmp


#----------------- Filter out reads which are associated with repeat element identified in the mouse genome --------------------#

# Take the result from Repeat Masker analysis on the annotation of the repeats that are present in the query sequence
## Tool: repeat masker
## Species: mouse genome
## Note: done before!
#sed -i '1,3d' 1/chr1.fa.out 2/chr2.fa.out 3/chr3.fa.out 4/chr4.fa.out 5/chr5.fa.out 6/chr6.fa.out 7/chr7.fa.out 8/chr8.fa.out 9/chr9.fa.out 10/chr10.fa.out 11/chr11.fa.out 12/chr12.fa.out 13/chr13.fa.out 14/chr14.fa.out 15/chr15.fa.out 16/chr16.fa.out 17/chr17.fa.out 18/chr18.fa.out 19/chr19.fa.out X/chrX.fa.out Y/chrY.fa.out
#cat 1/chr1.fa.out 2/chr2.fa.out 3/chr3.fa.out 4/chr4.fa.out 5/chr5.fa.out 6/chr6.fa.out 7/chr7.fa.out 8/chr8.fa.out 9/chr9.fa.out 10/chr10.fa.out 11/chr11.fa.out 12/chr12.fa.out 13/chr13.fa.out 14/chr14.fa.out 15/chr15.fa.out 16/chr16.fa.out 17/chr17.fa.out 18/chr18.fa.out 19/chr19.fa.out X/chrX.fa.out Y/chrY.fa.out > mm10_RepeatMasker.fa.out


# Prepare the BED file with only three columns: chr, start and end of each repeat regions that correspond to the repeat elements in the mouse genome
awk '{print $5 "\t" $6 "\t" $7}' $repeat_elements > $workdir/1.mapping_cleaning/output/mouse_repeat_regions_to_remove.bed

# Remove reads where overlap is found
## Tool: bedtools intersect
## Options:
### -v: keep only reads in the AAV mapping file that have no overlap in the regions identified in the BED file
### -abam: each read (BAM alignment) in the AAV mapping file is compared to the regions in the BED file in search of overlaps
### -b: to specify the following BED file
../bin/bedtools intersect -v -abam $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.bam -b $workdir/1.mapping_cleaning/output/mouse_repeat_regions_to_remove.bed > $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse.no_repeat_elements.bam 
