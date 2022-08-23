#!/bin/bash

#-------------------------- Set up the parameters -----------------------#

# Specify the working directory
workdir="$HOME/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/mark_duplicates"

# Path to the alignment file ("output.bam")
path_to_file="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/results"
path_to_destination="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/mark_duplicates/picard_mark_duplicates"



#-------------------------- Run the program -----------------------#

## Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"

    # Subset the marked pcr duplicates reads of each sample
    #samtools view -b --TAGGING_POLICY -f 1024 -@ 12 $path_to_file/${sample}_AAV_corrected/align_AAV/output.bam > $path_to_destination/${sample}_output.pcr_duplicates.bam
    java -jar ~/picard/build/libs/picard-2.26.2-1-g3a28787-SNAPSHOT-all.jar MarkDuplicates REMOVE_DUPLICATES=false TAGGING_POLICY=All \
      I=$path_to_file/${sample}_AAV_corrected/align_AAV/align_AAV.sorted.bam \
      O=$path_to_destination/${sample}_align_AAV.sorted.marked_duplicates.bam \
      M=$path_to_destination/${sample}_marked_dup_metrics.txt

done <$workdir/sample_list.txt
