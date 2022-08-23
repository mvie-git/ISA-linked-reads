#!/bin/bash

# Version 1.0
# This script was last modified on 19/11/2021
# Script Goal: Prepare IS data file from chimeric molecules for karyoploteR package on VECTOR genome
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: 
# Extract useful information on each IS from mouse molecule BED file after identify chimeric molecules
# 
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the mouse molecule BED file of interest in quotes
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
	mkdir -p $workdir/visualization/karyoploteR/input/vector/IS_chimeric_molecules/short

	# Give the file name of the AAV mapping BAM file filtered (no false positives)
	molecule_file="${workdir}/4.identify_chimeric_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.chimeric.bed"


	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################

	#----------------- Extract chr, start and end positions of each IS identified for karyoploteR package --------------------#

	# Extract chr, start and end positions for each alignment
	awk '{print $1 "\t" $2 "\t" $3}' $molecule_file > $workdir/visualization/karyoploteR/input/vector/IS_chimeric_molecules/short/${sample}_AAV_chimeric_IS.tmp

	# Rename the first column in 'chr'
	sed 1d $workdir/visualization/karyoploteR/input/vector/IS_chimeric_molecules/short/${sample}_AAV_chimeric_IS.tmp | sed '1ichr\tstart\tend' > $workdir/visualization/karyoploteR/input/vector/IS_chimeric_molecules/short/${sample}_AAV_chimeric_IS.bed

	# Clean
	rm $workdir/visualization/karyoploteR/input/vector/IS_chimeric_molecules/short/${sample}_AAV_chimeric_IS.tmp

	# Launch R script "karyoploteR_vector_chimeric_molecules_IS.R"

	#------------------------------------------- END --------------------------------------------------#

done < sample_list.txt