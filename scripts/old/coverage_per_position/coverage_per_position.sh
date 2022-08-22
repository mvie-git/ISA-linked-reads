#!/bin/bash

#-------------------------- Set up the parameters -----------------------#

# Specify the working directory
workdir="$HOME/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/coverage_per_position"

# Path to the alignment file ("output.bam")
path_to_file="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/results"
path_to_destination="/home/mallaury/Bureau/phd/tell_seq/PROJECT_GSDIa/tellseq_analysis/coverage_per_position"



#-------------------------- Run the program -----------------------#

# Create the output file to put the results from this script
# If doesn't exist
if [ ! -e "$path_to_destination/coverage_per_position_all_sample.txt" ] ; then
    echo -e 'sample,chromosome,position,coverage' > "$path_to_destination/coverage_per_position_all_sample.csv"
fi

## Looping through the content of a sample index list file: sample_list.txt
while read sample; do

  echo "<============ Starting sample n°$sample ============>"

    # Calculate coverage depth for all the position of the vector sequence (-a option include zero coverage positions)
    samtools depth -a $path_to_file/${sample}_AAV_corrected/align_AAV/output.bam > coverage_per_position_sample.tmp

    # Add which sample is processing in the first column & convert to csv format
    awk '{print '${sample}', $0}' coverage_per_position_sample.tmp | awk '{print $1 "," $2 "," $3 "," $4}' >> "$path_to_destination/coverage_per_position_all_sample.csv"

done <$workdir/sample_list.txt

# Clean up!
rm *.tmp