#!/bin/bash

#-------------------------- Set up the parameters -----------------------#

# Specify the working directory
workdir="$HOME/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/mark_duplicates"

# Path to the alignment file ("output.bam")
path_to_file="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/results"
path_to_destination="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/mark_duplicates/clumpify"



#-------------------------- Run the program -----------------------#

## Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"

    # Prepare input file: convert bam file to fastq
    bedtools bamtofastq -i $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.bam -fq $path_to_destination/${sample}_align_AAV.sorted.fq

    # Compress files
    pigz $path_to_destination/${sample}_align_AAV.sorted.fq

    # Remove exact duplicates
    clumpify.sh in=$path_to_destination/${sample}_align_AAV.sorted.fq.gz out=$path_to_destination/${sample}_align_AAV.sorted_clumped.fq.gz dedupe subs=0

    # # Mark exact duplicates but don't remove them (they get " duplicate" appended to the name)
    # clumpify.sh in=$path_to_destination/${sample}_align_AAV.sorted.fq.gz out=$path_to_destination/${sample}_align_AAV.sorted_clumped.fq.gz markduplicates subs=0

    # # Remove duplicates, allowing up to 5 substitutions between copies
    # clumpify.sh in=$path_to_destination/${sample}_align_AAV.sorted.fq.gz out=$path_to_destination/${sample}_align_AAV.sorted_clumped.fq.gz dedupe subs=5

    # # Remove optical duplicates only (duplicates within 40 pixels of each other - this value is platform-specific)
    # clumpify.sh in=r$path_to_destination/${sample}_align_AAV.sorted.fq.gz out=$path_to_destination/${sample}_align_AAV.sorted_clumped_optical.fq.gz dedupe optical dist=40 spantiles=f


done <$workdir/sample_list.txt


# Clean
