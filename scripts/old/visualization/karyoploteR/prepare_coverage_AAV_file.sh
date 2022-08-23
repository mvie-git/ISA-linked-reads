#!/bin/bash

# Version 1.0
# This script was last modified on 16/11/2021
# Script Goal: Prepare coverage file for karyoploteR package on VECTOR genome
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: 
# Calcul coverage data on the AAV mapping BAM after filtering out false positives
# 
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file of interest in quotes
#
# Output: Coverage data files for Rstudio script


#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the number of threads to use
nb_threads=10

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/visualization/karyoploteR/input/vector/AAV_ITRs/5_ITR

	# Give the file name of the AAV mapping BAM file filtered (no false positives)
	mapping_file="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/ITR_searching/output/${sample}/${sample}_5_ITR.sorted.bam"


	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	#----------------- Calculate coverage data file for karyoploteR package --------------------#

	# Tool: bedtools genomecov
	## Options:
	### -d: report the depth at each genome position
	### -ibam: input file in BAM format
	../../bin/bedtools genomecov -d -ibam $mapping_file > $workdir/visualization/karyoploteR/input/vector/AAV_ITRs/5_ITR/${sample}_vector_coverage_position.txt

	# Launch R script "karyoploteR_vector_coverage.R"


	#------------------------------------------- END --------------------------------------------------#

done < sample_list.txt