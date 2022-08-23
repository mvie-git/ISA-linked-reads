#!/bin/bash

#-------------------------- Set up the parameters -----------------------#

# Specify the working directory
workdir="$HOME/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/similarity_vector_mouse_genome/subset_bam"

# Path to the alignment file ("output.bam")
path_to_file="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/results"
path_to_destination="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/similarity_vector_mouse_genome/subset_bam"



#-------------------------- Run the program -----------------------#

# # Create the output file to put the results from this script
# # If doesn't exist
# if [ ! -e "$path_to_destination/summary.txt" ] ; then
#     echo -e 'sample\tnb_reads_before\tnb_reads_after' > "$path_to_destination/summary.txt"
# fi

## Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"

    # Subset the reads associated with similar regions in the mouse genome (for each sample)
    bedtools intersect -v -abam $path_to_file/${sample}_AAV_corrected/align_AAV/output.bam -b ../BLAST/blast_output_vector_filter.bed > $path_to_file/${sample}_AAV_corrected/align_AAV/output.filter.similar_mouse.bam

    # # Doesn't work!
    # # Number of reads before filtering
    # nb_before_filter=samtools view $path_to_file/${sample}_AAV_corrected/align_AAV/output.bam | wc -l

    # # Number of reads after filtering
    # nb_after_filter=samtools view $path_to_file/${sample}_AAV_corrected/align_AAV/output.filter.similar_mouse.bam | wc -l

    # # Number of reads removed after filtering
    # nb_reads_removed=$((nb_before_filter-nb_after_filter))

    # # Fill the summary.txt file
    # echo -e "${sample}\t$nb_before_filter\t$nb_reads_removed" >> "$path_to_destination/summary.txt"


done <$workdir/sample_list.txt