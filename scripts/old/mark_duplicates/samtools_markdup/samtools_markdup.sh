#!/bin/bash

#-------------------------- Set up the parameters -----------------------#

# Specify the working directory
workdir="$HOME/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/mark_duplicates"

# Path to the alignment file ("output.bam")
path_to_file="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/results"
path_to_destination="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/mark_duplicates/samtools_markdup"
nb_threads=100



#-------------------------- Run the program -----------------------#

## Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"

    # The first sort can be omitted if the file is already name ordered
    samtools sort -n -o $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.tmp $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.bam

    # Add ms and MC tags for markdup to use later
    samtools fixmate -m $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.tmp $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.tmp

    # Markdup needs position order
    samtools sort -o $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.positionsort.tmp $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.tmp

    # Finally mark duplicates
    samtools markdup $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.positionsort.tmp $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.positionsort.markdup.bam


done <$workdir/sample_list.txt


# Clean
rm $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.tmp 
rm $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.tmp
rm $path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.name.fixmate.positionsort.tmp

