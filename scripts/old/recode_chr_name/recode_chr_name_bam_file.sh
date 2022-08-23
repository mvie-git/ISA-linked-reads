#!/bin/bash

# Version 1.0
# This script was last modified on 28/10/2021
# Script Goal: Recode chromosome name in BAM alignment files
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: Before Tell-Read analysis, the genome is indexed by the script 'generateGenomeIndexBed' but it sort the chromosome names by the size of it and not by the conventional names
# Extract header and body of the BAM alignment file to convert it to SAM format
# Apply a python script to recode all the chromosome names
# Concatenate the header and the modified body SAM files
# Convert to BAM format
# Sort BAM file (because the order of chromosomes have changed but it's the same than in the header BAM file!)
# Index BAM file
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this BAM alignment file is stored
# 2. Write the name of the AAV/mouse mapping BAM file of interest in quotes
#
# Output: The script generates a modified BAM file with the right chromosome names
# 
# Concept:
# Recode all chromosome names of a BAM file due to an error during the Tell-Seq analysis
# 
# Note:
# Files from the mouse mapping are heavy! Check the time (it may be long)
# 
# Source: https://www.biostars.org/p/46104/



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/recode_chr_name/)
workdir="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Create the output directory (only if it doesn't exist already!)
mkdir -p $workdir/recode_chr_name/output

# Set the number of threads to use
nb_threads=10

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

    echo "*** Start $sample ***"

	# Give the file name of the mouse mapping BAM file of interest
	mapping_file="/media/mallaury/WLG/151/align_mouse/align_mouse_temp/align_mouse.sorted.bam"
	index_file="/media/mallaury/WLG/151/align_mouse/align_mouse_temp/align_mouse.sorted.bam"

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/recode_chr_name/output/${sample}

	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	# INFO! If you doesn't want to decompress the BAM into SAM file, you need to extract the header of the samtools view command (SAM file) to put it back before convert it back to BAM format


	#----------------- Decompress the header of the BAM file into SAM file --------------------#

	# Save the header lines 
	../bin/samtools view -H $mapping_file > $workdir/recode_chr_name/output/${sample}/${sample}_header.sam

	# Execute the python script to recode all chromosome names
	## Argument: header SAM file is the input file to be modified
	python3 recode_chr_name.py $workdir/recode_chr_name/output/${sample}/${sample}_header.sam



	#-------------- Replace the header in BAM file with the modified header ---------------#

	# Combine header and body modified
	../bin/samtools reheader $workdir/recode_chr_name/output/${sample}/${sample}_header.sam $mapping_file > $workdir/recode_chr_name/output/${sample}/${sample}_chr_recode.bam

	# Index the BAM file 
	../bin/samtools index $workdir/recode_chr_name/output/${sample}/${sample}_chr_recode.bam


	# Clean up!
	## Replace the BAM file
	#mv $workdir/recode_chr_name/output/${sample}/${sample}_chr_recode.bam $mapping_file

	## Replace the index file
	#mv $workdir/recode_chr_name/output/${sample}/${sample}_chr_recode.bam.bai $index_file


	#------------------------------------------- END --------------------------------------------------# 

done < sample_list.txt