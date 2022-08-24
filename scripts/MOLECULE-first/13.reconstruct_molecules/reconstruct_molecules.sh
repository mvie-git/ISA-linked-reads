#!/bin/bash

# Version 1.0
# This script was last modified on 19/07/2022 (UPDATE)
# Script Goal: Get molecules
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After filtering reads from mapping results (to both genomes), identify "molecules" recovered
# From the mapping BAM file filtered with only mapped reads
# Recover molecules with the script from 10X Genomics (+ adapted for Tell-Seq analysis)
# Get information on molecules: chr, start, end, barcode sequence, length
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file (filtered) of interest in quotes
# 3. Write the name of the script used to recover molecules from 10X Genomics
#
# Output: The script generates two BED molecule files
# 1. A BED file with information on molecules recovered from mapped and filtered reads: chr, start, end, barcode sequence, length
#
# Concept:
# Create a BED molecule file with all the necessary informations from AAV/mouse mapping BAM file filtered



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder
workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping" #----------------------> !CHANGE IT!#

# Set the name of the sample
#sample=149 #------------------------> !CHANGE IT!#

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/get_molecules/${sample}

	# Give the file name of the mapping BAM file
	input_BAM=${workdir}/add_barcode_to_BAM/11.add_barcode_to_BAM_file/${sample}/${sample}_subset_BAM_filtered.BX_tag.bam
	
	# Give the file name of the script used to get the molecules (10X Genomics & Tell-Seq)
	get_molecules_script="${workdir}/get_molecules/bin/GetMoleculesInfo_10X_modified.py"

	# Set the number of threads to use
	nb_threads=12


	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################


	#---------------------------------- Prepare the input BAM file to keep only mapped reads --------------------------------#

	# # Extract the header of the input BAM file
	# samtools view -H $input_BAM > ${sample}_header.sam

	# # Extract mapped reads from the input BAM file
	# samtools view -F 4 $input_BAM > ${sample}_body.sam

	# # Concatenate both sam files
	# cat ${sample}_header.sam ${sample}_body.sam > ${sample}_mapped.sam

	# # Convert SAM to BAM format
	# samtools view -Sb ${sample}_mapped.sam > ${sample}_mapped.bam

	# # Sort and index the new BAM file
	# samtools sort -@ $nb_threads ${sample}_mapped.bam > ${sample}_mapped.sorted.bam
	# samtools index ${sample}_mapped.sorted.bam

	# # Clean
	# rm *.sam ${sample}_mapped.bam ${sample}_mapped.sorted.bam ${sample}_mapped.sorted.bam.bai


	#---------------------------------- Create and modify the "molecule.tsv" file from AAV mapping BAM file filtered --------------------------------#


	# (If it's not already done!) Index the BAM alignment files before use the get molecules script
	#samtools index -@ $nb_threads $input_BAM 

	# Create the corresponding "molecule.tsv" file from the BAM file
	python3 $get_molecules_script $input_BAM ${sample}.molecules.tmp

	# Add a column for the length of each molecule
	## Options: 
	### Set tab delimiter 'OFS'
	awk 'BEGIN {OFS="\t"} {$6=$3-$2} 1' ${sample}.molecules.tmp > ${sample}.molecules_2.tmp

	# Add a header: chromosome, start, end, barcode sequence, number of reads, length of the molecule
	sed -e '1i\chromosome\tstart\tend\tbarcode\tnb_reads\tlength' ${sample}.molecules_2.tmp > $workdir/get_molecules/${sample}/${sample}.molecules.bed

	# Clean
	rm *.molecules.tmp *.molecules_2.tmp


	#----------------------------------------------------------------- END ----------------------------------------------------------------#

done < sample_list.txt
