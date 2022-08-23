#!/bin/bash

# Version 1.0
# This script was last modified on 18/11/2021
# Script Goal: Prepare IS data file from chimeric reads for karyoploteR package on MOUSE genome
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: 
# Extract useful information on each IS from mouse mapping BAM file after identify chimeric reads
# 
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the mouse mapping BAM file of interest in quotes
#
# Output: IS data files for Rstudio script


#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/home/mallaury/Bureau/phd/PROJECT/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the number of threads to use
nb_threads=10

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/visualization/karyoploteR/input/mouse/IS_chimeric_reads/short

	# Give the file name of the AAV mapping BAM file filtered (no false positives)
	mapping_file="${workdir}/2.identify_chimeric_reads/output/short/${sample}/${sample}_mouse_mapping_filtered.chimeric.bam"


	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	#----------------- Extract chr, start and end positions of each IS identified for karyoploteR package --------------------#

	# Extract length, chr and end positions for each alignment
	../../bin/samtools view -@ $nb_threads $mapping_file | awk '{print $2 "\t" $3 "\t" $4}' > $workdir/visualization/karyoploteR/input/mouse/IS_chimeric_reads/short/file1.tmp

	# Change header
	sed 1d $workdir/visualization/karyoploteR/input/mouse/IS_chimeric_reads/short/file1.tmp | sed '1ilength\tchr\tstart' > $workdir/visualization/karyoploteR/input/mouse/IS_chimeric_reads/short/${sample}_mouse_chimeric_IS.bed
	
	# Launch R script "karyoploteR_vector_chimeric_reads_IS.R"

	# Clean!
	rm $workdir/visualization/karyoploteR/input/mouse/IS_chimeric_reads/short/file1.tmp


	#------------------------------------------- END --------------------------------------------------#

done < sample_list.txt