#!/bin/bash

# Version 1.0
# This script was last modified on 19/10/2021
# Script Goal: Chimeric reads
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After filtering the output from AAV and mouse mapping results, identify which reads are "chimeric" (meaning mapping to both genomes AAV/mouse)
# Extract list of unique reads present in each genome (AAV/mouse)
# Compare the two lists in order to identify "chimeric" reads
# Subset the AAV/mouse mapping BAM file (filtered) in order to keep only information on chimeric reads
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file (filtered) of interest in quotes
# 3. Write the name of the mouse mapping BAM file (filtered) of interest in quotes
#
# Output: The script generates two BAM files
# 1. List of unique reads present in both genomes (AAV/mouse) = "chimeric" reads
# 2. Subset of the AAV mapping BAM file (filtered) with only "chimeric" reads
# 3. Subset of the mouse mapping BAM file (filtered) with only "chimeric" reads
# 
# Concept:
# Extract unique reads from AAV mapping BAM file (filtered) and compare it to the list of unique reads from mouse mapping BAM file (filtered)
# 
# Note:
# Files from the mouse mapping are heavy!



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/2.identify_chimeric_reads/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
#sample=149 #------------------------> !CHANGE IT!#

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create the output directory (only if it doesn't exist already!)
	mkdir -p $workdir/2.identify_chimeric_reads/output/short

	# Give the file name of the AAV/mouse mapping BAM file (filtered)
	AAV_mapping_file="${workdir}/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam"
	mouse_mapping_file="${workdir}/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam"

	# Set the number of threads to use
	nb_threads=10

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/2.identify_chimeric_reads/output/short/${sample}

	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	#----------------- Extract list of unique reads present in both genomes (AAV/mouse) = "chimeric" reads --------------------#

	# Extract list of unique reads mapped to AAV after filtering false positives
	../bin/samtools view -@ nb_threads $AAV_mapping_file | awk '{print $10}' | sort | uniq > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_reads_after_filtering.txt
	wc -l $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_reads_after_filtering.txt

	# Extract list of unique reads mapped to mouse after filtering false positives
	../bin/samtools view -@ nb_threads $mouse_mapping_file | awk '{print $10}' | sort | uniq > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_mouse_reads_after_filtering.txt
	wc -l $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_mouse_reads_after_filtering.txt


	# Extract reads present in only AAV mapping results and reads common in both genomes (-2: supress column with reads only in mouse genome)
	comm -2 $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_reads_after_filtering.txt $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_mouse_reads_after_filtering.txt > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_after_filtering_comm_output.tmp

	# Extract the first column with reads only present in AAV mapping BAM file filtered (meaning AAV only reads)
	awk -F "\t" '{print $1}' $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_after_filtering_comm_output.tmp > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_AAV_only.txt

	# Remove empty lines and count the number of reads mapped to only AAV genome (meaning AAV only reads)
	sed -i '/^$/d' $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_AAV_only.txt 
	wc -l $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_AAV_only.txt 

	# Extract the second column with reads present in both genomes (AAV and mouse)
	awk -F "\t" '{print $2}' $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_after_filtering_comm_output.tmp > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_chimeric.txt

	# Remove empty lines and count the number of reads mapped to both genomes (AAV and mouse)
	sed -i '/^$/d' $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_chimeric.txt
	wc -l $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_chimeric.txt

	# Clean up!
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_after_filtering_comm_output.tmp
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_reads_after_filtering.txt # Takes a lot of disk place
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_mouse_reads_after_filtering.txt # Takes a lot of disk place



	#-------- Subset the AAV mapping BAM file (filtered) in order to keep only information on chimeric reads --------#

	# INFO! If you doesn't want to decompress the BAM into SAM file, you need to extract the header of the samtools view command (SAM file) to put it back before convert it back to BAM format

	# Save the header lines 
	../bin/samtools view -H $AAV_mapping_file > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_header_AAV_mapping_filtered.sam

	# Extract mapped reads from the BAM filtered file on the AAV genome with the list of "chimeric" reads
	## Filter alignments using list of "chimeric" reads (-w option for exact string match)
	../bin/samtools view -@ $nb_threads $AAV_mapping_file | grep -w -f $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_chimeric.txt > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_body_AAV_mapping_filtered.sam

	# Combine header and body 
	cat $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_header_AAV_mapping_filtered.sam $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_body_AAV_mapping_filtered.sam > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.sam 

	# Convert output.filtered.sam to BAM format 
	../bin/samtools view -@ $nb_threads -b $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.sam  > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.bam

	# Sort and index the BAM file 
	../bin/samtools sort -@ $nb_threads $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.bam -o $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_AAV_mapping_filtered.chimeric.bam
	../bin/samtools index $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_AAV_mapping_filtered.chimeric.bam

	# Clean up!
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_header_AAV_mapping_filtered.sam
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_body_AAV_mapping_filtered.sam
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.sam
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_AAV_mapping_filtered.bam



	#-------- Subset the mouse mapping BAM file (filtered) in order to keep only information on chimeric reads --------#

	# INFO! If you doesn't want to decompress the BAM into SAM file, you need to extract the header of the samtools view command (SAM file) to put it back before convert it back to BAM format

	# Save the header lines 
	../bin/samtools view -H $mouse_mapping_file > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_header_mouse_mapping_filtered.sam

	# Extract mapped reads from the BAM filtered file on the mouse genome with the list of "chimeric" reads
	## Filter alignments using list of "chimeric" reads (-w option for exact string match)
	../bin/samtools view -@ $nb_threads $mouse_mapping_file | grep -w -f $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_list_AAV_mouse_reads_chimeric.txt > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_body_mouse_mapping_filtered.sam

	# Combine header and body 
	cat $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_header_mouse_mapping_filtered.sam $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_body_mouse_mapping_filtered.sam > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.sam 

	# Convert output.filtered.sam to BAM format 
	../bin/samtools view -@ $nb_threads -b $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.sam  > $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.bam

	# Sort and index the BAM file 
	../bin/samtools sort -@ $nb_threads $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.bam -o $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam
	../bin/samtools index $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam

	# Clean up!
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_header_mouse_mapping_filtered.sam
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_body_mouse_mapping_filtered.sam
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.sam
	rm $workdir/2.identify_chimeric_reads/output/short/${sample}/${sample}_output_all_mouse_mapping_filtered.bam


	#------------------------------------------- END --------------------------------------------------# 

done < sample_list.txt
