#!/bin/bash

# Version 1.0
# This script was last modified on 03/11/2021 (UPDATE? )
# Script Goal: Get molecules
# Author: Mallaury Vie
# Contributions by: 
# 
# Goal: After filtering reads from mapping results (to both genomes), identify "molecules" recovered
# From the AAV mapping BAM file filtered
# Recover molecules with the script from 10X Genomics (+ adapted for Tell-Seq analysis)
# Get information on molecules: chr, start, end, barcode sequence, length
#
# Input: Specify the following parameters
# 1. Set the path to the directory where this file is stored
# 2. Write the name of the AAV mapping BAM file (filtered) of interest in quotes
# 3. Write the name of the mouse mapping BAM file (filtered) of interest in quotes
# 4. Write the name of the script used to recover molecules from 10X Genomics
#
# Output: The script generates two BED molecule files
# 1. A BED file with information on molecules recovered from AAV mapped and filtered reads: chr, start, end, barcode sequence, length
# 2. A BED file with information on molecules recovered from mouse mapped and filtered reads: chr, start, end, barcode sequence, length
#
# Concept:
# Create a BED molecule file with all the necessary informations from AAV/mouse mapping BAM file filtered
# 
# Note:
# Files from the mouse mapping are heavy!



#----------------- CHANGES: Set parameters in this section manually --------------------#

# Set the directory of the script as the working folder (e.g /home/user/studyname/ISA/3.get_molecules/)
workdir="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA" #----------------------> !CHANGE IT!#

# Set the name of the sample
#sample=149 #------------------------> !CHANGE IT!#

#ALTERNATIVE: Looping through the content of a sample index list file: sample_list.txt #----------------> !Change IT!#
while read sample; do

	echo "*** Start $sample ***"

	# Create the output directory (only if it doesn't exist already!)
	mkdir -p $workdir/3.get_molecules/output/short

	# Create a folder to store all the results for this sample
	mkdir -p $workdir/3.get_molecules/output/short/${sample}

	# Give the file name of the AAV/mouse mapping BAM file (filtered)
	AAV_mapping_file="${workdir}/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam"
	mouse_mapping_file="${workdir}/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam"

	# Give the file name of the script used to get the molecules (10X Genomics & Tell-Seq)
	get_molecules_script="${workdir}/3.get_molecules/input/GetMoleculesInfo_10X_modified.py"

	# Set the number of threads to use
	nb_threads=10


	#----------------- END: No changes are needed below this line --------------------#



	##################################################################################
	######                             Main Script                              ######
	##################################################################################


	#---------------------------------- Create and modify the "molecule.tsv" file from AAV mapping BAM file filtered --------------------------------#


	# Index the BAM alignment files before use the get molecules script
	../bin/samtools index -@ $nb_threads $AAV_mapping_file $workdir/1.mapping_cleaning/output/vector/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam.bai

	# Create the corresponding "molecule.tsv" file from the AAV mapping BAM file (filtered)
	python3 $get_molecules_script $AAV_mapping_file $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.tmp

	# Add a column for the length of each molecule
	## Options: 
	### Set tab delimiter 'OFS'
	awk 'BEGIN {OFS="\t"} {$6=$3-$2} 1' $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.tmp > $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.tmp2

	# Add a header: chromosome, start, end, barcode sequence, number of reads, length of the molecule
	sed -e '1i\chromosome\tstart\tend\tbarcode\tnb_reads\tlength' $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.tmp2 > $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.bed

	# Clean
	rm $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_AAV.molecules.tmp*


	#---------------------------------- Create and modify the "molecule.tsv" file from mouse mapping BAM file filtered --------------------------------#


	# Index the BAM alignment files before use the get molecules script
	../bin/samtools index -@ $nb_threads $mouse_mapping_file $workdir/1.mapping_cleaning/output/mouse/${sample}/${sample}_mapping.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.no_similar_mouse_short.no_repeat_elements.bam.bai

	# Create the corresponding "molecule.tsv" file from the mouse mapping BAM file (filtered)
	python3 $get_molecules_script $mouse_mapping_file $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.tmp

	# Add a column for the length of each molecule
	## Options: 
	### Set tab delimiter 'OFS'
	awk 'BEGIN {OFS="\t"} {$6=$3-$2} 1' $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.tmp > $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.tmp2

	# Add a header: chromosome, start, end, barcode sequence, number of reads, length of the molecule
	sed -e '1i\chromosome\tstart\tend\tbarcode\tnb_reads\tlength' $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.tmp2 > $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.bed

	# Compress BED file
	pigz $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.bed

	# Clean
	rm $workdir/3.get_molecules/output/short/${sample}/${sample}_mapping_mouse.molecules.tmp*


	#----------------------------------------------------------------- END ----------------------------------------------------------------#

done < sample_list.txt
