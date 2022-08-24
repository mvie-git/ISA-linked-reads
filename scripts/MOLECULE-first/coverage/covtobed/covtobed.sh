# Workdir
workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping"

# BAM file
bam_file="$workdir/output/152_mapping_mouseANDvector.sorted.bam"

# Output file name
depth_file="$workdir/covtobed/152_mapping_mouseANDvector.depth"
genomecov_file="$workdir/covtobed/152_mapping_mouseANDvector.genomecov"
covtobed_file="$workdir/covtobed/152_mapping_mouseANDvector.covtobed"


#------------ Extracting coverage information with samtools -------------#

## Samtools depth is a tool that calculates the sequence coverage at each position
# Output format: chr, position, coverage
#samtools depth -a $bam_file | pigz > ${depth_file}.txt.gz

## Bedtools genomecov is a tool that produce a standard BED file (chr, start, end, coverage)
# Adjacent bases with the same coverage will be collapsed in the same interval
bedtools genomecov -ibam $bam_file -bga | pigz > ${genomecov_file}.bed.gz


#--------- Getting the regions covered within a coverage range ------------#

## Covtobed is a tool to produce a coverage profile, that takes a sorted BAM file as input
## and will produce a BED file. Optionnally specifying to only prints those intervals 
## having a coverage within a range (or only the regions with no coverage)

# Install covtobed
conda install -y -c bioconda covtobed

# Only outputs the regions with no coverage
covtobed -m 0 -M 1 $bam_file | pigz > ${covtobed_file}.not_covered.bed.gz

# Only outputs the regions above a specific coverage (example: the mean/median)


#------------ Some statistics on the BAM file -----------#

## Calculate average genome coverage on aligned BAM files
# Be sure to include regions with no coverage
zcat ${depth_file}.txt.gz | awk '{sum+=$3} END { print "Average = ",sum/NR}'

