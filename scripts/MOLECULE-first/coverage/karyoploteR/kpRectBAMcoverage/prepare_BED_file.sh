#!/bin/bash

# Version 1.0
# This script was last modified on 20/07/2021
# Script Goal: Prepare input data file from BAMfile for karyoploteR package
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: 
# Extract useful information on each read from a mapping BAM file
# 
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the mapping BAM file of interest in quotes
#
# Output: BED input data file for Rstudio script


#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/visualization/)
workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/karyoploteR/kpRectBAMcoverage" #----------------------> !CHANGE IT!#

sambamba=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/bin/sambamba

# Set the number of threads to use
nb_threads=12

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/input

	# Give the file name of the mapping BAM file
	input_BAM=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/${sample}/${sample}_mapping_mouseANDvector.sorted.chrVirus.mapped.bam


	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################


	# Identify and remove duplicate reads (PCR and optical duplicates) from the AAV mapping BAM file
	## Tool: sambamba markdup
	## Options:
	### -r: remove duplicates instead of just marking them
	### -t: number of threads
	### -p: show progressbar in the terminal
	$sambamba markdup -r -t $nb_threads -p $input_BAM /media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/${sample}/${sample}_mapping_mouseANDvector.sorted.chrVirus.mapped.no_duplicates.bam



	#----------------- Extract chr, start and end positions of each read in the BAM file for karyoploteR package --------------------#

	# # Extract length, chr and start positions for each alignment
	# samtools view -@ $nb_threads $input_BAM | awk '{print $2 "\t" $3 "\t" $4}' > ${workdir}/input/file1.tmp

	## Alternative: bedtools bamtobed -i reads.bam | head -3
	bedtools bamtobed -i /media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/${sample}/${sample}_mapping_mouseANDvector.sorted.chrVirus.mapped.no_duplicates.bam > ${workdir}/input/${sample}.tmp
	#bedtools bamtobed -i $input_BAM > ${workdir}/input/${sample}.tmp

	# Calcul length
	awk '{print $1 "\t" $2 "\t" $3}' ${workdir}/input/${sample}.tmp | awk '$4=$3-$2 {print $0}' | awk '{print $4 "\t" $1 "\t" $2}' > ${workdir}/input/${sample}.tmp2

	# Change header
	sed 1d ${workdir}/input/${sample}.tmp2 | sed '1ilength\tchr\tstart' > ${workdir}/input/${sample}.bed
	
	# Launch R script "karyoploteR.R"

	# Clean!
	rm ${workdir}/input/${sample}.tmp ${workdir}/input/${sample}.tmp2


	#------------------------------------------- END --------------------------------------------------#

done < sample_list.txt