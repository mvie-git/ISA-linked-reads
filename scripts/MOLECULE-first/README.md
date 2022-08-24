**Molecule-first** strategy consists in gathering all reads sharing the same barcode into molecules then identify chimeric molecules.
There is three steps:
1. Gather all reads sharing the same barcode to reconstruct the initial DNA molecule
2. Map all the molecules against the hybrid genome
3. Identify chimeric molecules meaning with reads mapped to both genomes


## Generate TELL-Seq linked-reads from sequencing data

> The Tell-Read module of the [TELL-Seq software](https://www.universalsequencing.com/) demultiplex and correct error barcoded reads from the initial raw sequencing data.

The Tell-Seq pipeline (`Tell-Read` module) takes as input the initial raw sequencing FASTQ data: 
- **demultiplexing**: they are further demultiplexed into sample separated I1 (index 1) reads, R1 reads, R2 reads, based on I2 (index 2) reads. I1 reads are the TELL_seq barcode sequences. For each sequencing library construction, a set of unique barcode sequences was randomly chosen from a 2.4 billion-barcode pool. These sample-demultiplexed FASTQ files are saved as the raw data output files. 
- **adaptor trimming**: adapter sequences are then trimmed using *cutadapt* utility.
- **error barcoded reads correction**: the unique barcodes associated with only one read are most likely caused by sequencing errors in the barcode. These barcodes are first identified if they are 1-base mismatches with one of the barcode associated with multiple reads, and then error-corrected. Barcodes with errors after this step are filtered out. The erroneous barcodes along with their associated reads are removed and excluded from the rest of analyses. The remainin R1 and R2 reads, along with their associated I1 reads (barcodes) are the TELL-seq linked reads. They are the input for downstream analysis

To run Tell-Read module:
```
# Variables
sample=149 # Name of the sample used
output_dir=demultiplex_vector # Name the folder where the result will be stored
genome_dir=${workdir}/genomes/vector # Path to the genome reference folder
genome_fasta=vector.fasta # Name of the fasta file of the reference genome
```
1. Prepare Docker image:
```
# Load the Docker image and check
docker load -i docker-tellread
docker images
```
2. Prepare genome directory
```
# Create directory for genome reference fasta files
mkdir genomes
cd genomes/
mkdir vector mouse

# This script generates the genome's indexes and BED files for each reference
$path/tellread-release/generateGenomeIndexBed.sh genomes/vector/vector.fasta
$path/tellread-release/generateGenomeIndexBed.sh genomes/mouse/mouse.fasta
```
3. Prepare raw sequencing data files
```
# (Optional) Check if the name of files correspond to the right content
## index 2 (8-base) | Convert I1 to I2 files for each lane
zcat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz | head
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I2.fastq.gz # first lane
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I2.fastq.gz # second lane

## index 1 (18-base) | Convert R2 to I1 files for each lane
zcat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz | head
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz # first lane
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I1.fastq.gz # second lane

## read 2 | Convert R3 to R2 files for each lane 
zcat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R3.fastq.gz | head
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R3.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz # first lane
mv $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R3.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R2.fastq.gz # second lane

# Merge lane 1 and lane 2 into a single file
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I2.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I2.fastq.gz
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.I1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.I1.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I1.fastq.gz
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R2.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R2.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R2.fastq.gz
cat $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.R1.fastq.gz $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.R1.fastq.gz > $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R1.fastq.gz

# Clean single lane files
rm $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane1.gcap_dev.*
rm $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane2.gcap_dev.*
```
4. Run Tell-Read script on FASTQ sequencing data files
```
$workdir/tellread-release/run_tellread_fq.sh \
-i1 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I1.fastq.gz \
-i2 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.I2.fastq.gz \
-r1 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R1.fastq.gz \
-r2 $workdir/${sample}/GC107${sample}.210106_WLG_GR002.210219.NovaSeq1.FCB.lane12.gcap_dev.R2.fastq.gz \
-o $workdir/${sample}/${output_dir} \
-f $workdir/${genome_dir}/ \
-s ${sample} \
-g ${genome_fasta}
```

To see the distribution of length reads:
```
zcat demultiplex_AAV_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'
```

## Prepare an index file of the TELL-seq linked reads
> It is useful to have the list of all the reads (1. read ID) and their corresponding barcode (2. barcode sequence) in a text file for downstream analysis. The script is called ![prepare_index_file.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/prepare_index_file)



## BWA: Map TELL-seq linked reads against the hybrid genome
> BWA is a software package for mapping low-divergent sequences against a large reference genome. BWA-MEM is better for longer sequences ranged from 70bp to 1Mbp.

![First step: Hybrid genome mapping](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/images/1_mapping_hybrid_genome.pdf "Hybrid genome mapping (MOLECULE-FIRST)")

**BWA-mem** takes as input the demultiplexed and error barcoded corrected TELL-seq linked reads (FASTQ) for mapping those reads against the hybrid genome. The two FASTA reference genome files are required to create an index of both references together with `BWA index` tool.
```
# Concatenate the two FASTA reference genome files
cat mouse.fasta vector.fasta > mouseANDvector.fasta

# Create index files for BWA with bwa index
${path}/bwa index -a bwtsw mouseANDvector.fasta
```
Then run `BWA mem`:
```
# Variables
nb_threads=1 # Number of threads to use
reference_genome=${path}/mouseANDvector.fasta # Hybrid genome reference FASTA file
R1_fastq_file=${path}/demultiplex_AAV_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz # Input FASTQ data file
R2_fastq_file=${path}/demultiplex_AAV_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz
output_dir=${path}/results/${sample} # Name the folder where the result will be stored
output_name=${sample}_mapping.bam

# Run bwa mem command | Sort and convert the output to BAM format
${path}/bwa mem -t ${nb_threads} $reference_genome $R1_fastq_file $R2_fastq_file | samtools view -bS - > ${output_name}
```

To view the distribution of the mapping reads against the hybrid genome with **IGV viewer**, it is required to sort the mapping BAM file by coordinates and index the output file:
```
# Sort the BAM file by coordinates (default) | Index the new BAM file
samtools sort -@ ${nb_threads} mapping.bam > mapping.sorted.bam
samtools index mapping.sorted.bam
```
Other useful commands to have some statistics on the alignment file with `samtools flagstat` command:
```
samtools flagstat mapping.sorted.bam
```

## First reads selection based on vector mapping results
![Second step: Vector selection mapping](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/images/2_vector_selection_.pdf "First reads selection based on vector (MOLECULE-FIRST)")

1. Clean the BAM mapping file to filter out unmapped reads, duplicates reads, multi-mapped reads and low quality mapped reads (![filter_input_BAM_file.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/1.filter_input_BAM_file_AAV))
2. Extract read ID of vector mapped reads (![select_readID_vector_selection.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/2.select_readID_vector_selection))
3. Subset the initial index file to retrieve barcodes of the list of read ID corresponding to vector selected mapped reads (![retrieve_barcodes_from_readID_vector.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/3.retrieve_barcodes_from_readID))
4. Retrieve all the reads IDs corresponding to the subsetting new index meaning retrieving all the reads sharing the same barcode as reads mapped to the vector sequence (![retrieve_all_readID_from_barcodes_vector.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/4.retrieve_all_readID_from_barcodes))
5. Subset the input BAM file to keep only reads matching the read ID in the new index (![subset_BAM_file_with_readID.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/5.subset_BAM_file_with_readID))


## Second reads selection based on mouse mapping results
![Third step: Mouse selection mapping](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/images/3_mouse_selection.pdf "Third reads selection based on mouse (MOLECULE-FIRST)")

6. Clean the BAM mapping file to filter out unmapped reads (mouse mapped reads only), duplicates reads, multi-mapped reads and low quality mapped reads. Here, we want to keep only « chimeric » molecules meaning reads mapped either on vector or mouse but sharing the same barcode. We filtered out barcodes associated with vector mapped reads only (![filter_input_BAM_file_mouse.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/6.filter_input_BAM_file_mouse))
7. Extract read ID of mouse mapped reads (![select_readID_mouse_selection.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/7.select_readID_of_interest))
8. Subset the initial index file to retrieve barcodes of the list of read ID corresponding to mouse selected mapped reads (![retrieve_barcodes_from_readID_mouse.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/8.retrieve_barcodes_from_readID))
9. Retrieve all the reads IDs corresponding to the subsetting new index meaning retrieving all the reads sharing the same barcode as reads mapped to the mouse genome (![retrieve_all_readID_from_barcodes_mouse.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/9.retrieve_all_readID_from_barcodes))
10. Subset the input BAM file to keep only reads matching the read ID in the new index (![subset_BAM_file_with_readID.sh](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/MOLECULE-first/10.subset_BAM_file_with_readID))


## Reconstruct chimeric molecules
![Fourth step: Reconstruct chimeric molecules](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/images/4_reconstruct_molecules.pdf "Reconstruct chimeric molecules (MOLECULE-FIRST)")

11. Add the barcode information as a tag in the new output BAM file with a tool named pysam which allow to parse a BAM file ![add_barcode_to_BAM_file.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/11.add_barcode_to_BAM_file/add_barcode_to_BAM_file.sh)
  > Info! For large BAM files (![add_barcode_to_large_BAM_file.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/11.add_barcode_to_BAM_file/with_large_BAM_file/add_barcode_to_large_BAM_file.sh))
12. Sort and index the BAM file (![sort_and_index_BAM_file.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/sort_and_index_BAM_file.sh))
13. Reconstruct molecules from a BAM file (![reconstruct_molecules.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/13.reconstruct_molecules/reconstruct_molecules.sh)) with adaptation of an original script from 10X Genomics (![GetMoleculesInfo_10X_modified.py](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/13.reconstruct_molecules/bin/GetMoleculesInfo_10X_modified.py))

Some useful commands to get some statistics on chimeric molecules:
```
# Average size of molecules
awk '{sum+=$5} END { print "Average = ",sum/NR}' chimeric.bed

# Number of chimeric molecules by chromosome
awk '{print $1}' chimeric.bed | sort | uniq -c
```
or some visualizations with **ggplot R package**...
- Boxplot: length of chimeric molecules, number of molecules per barcode
  - Preparation of the input BED file for ggplot R package (![prepare_data.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/ggplot/r_ggplot_boxplot/prepare_data.sh))
  - R script for boxplot function (![r_plot_boxplot.R](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/ggplot/r_ggplot_boxplot/r_plot_boxplot.R))
- Histogram: number of reads per molecule
  - Preparation of the input BED file for ggplot R package (![prepare_data.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/ggplot/r_ggplot_histogram/prepare_data.sh))
  - R script for histogram function (![r_plot_histo.R](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/ggplot/r_ggplot_histogram/r_plot_histo.R))


## Coverage

### Visualization with KaryoploteR package

- Distribution of vector mapped reads across the whole vector genome
  - Preparation of the BED input file for karyoplote R package (![prepare_BED_file.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/coverage/karyoploteR/kpRectBAMcoverage/prepare_BED_file.sh))
  - R script for kpRectBAMcoverage function (![kpRectBAMcoverage.R](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/coverage/karyoploteR/kpRectBAMcoverage/kpRectBAMcoverage.R))

- Distribution of the mouse mapped reads across the mouse genome
  - Preparation of the BED input file for karyoplote R function (![prepare_BED_file.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/coverage/karyoploteR/kpPlotDensity/prepare_BED_file.sh))
  - R script for kpPlotDensity function (![karyoploteR_density_mouse.R](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/coverage/karyoploteR/kpPlotDensity/karyoploteR_density_mouse.R))


### Visualization with PlotCoverage
> To see the profile of the coverage along the mouse genome: ![plotCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/plotCoverage.html) tool.

It is required to filter out unmapped reads, vector mapped reads and duplicates. The number of vector mapped reads is really low so it has a minor impact on a global overview of genome coverage from these data. But we need to filter out unmapped reads: some of them have coordinates because samtools will give them the same coordinated if their mate pair mapped somewhere on the reference genome.

```
samtools view -F 4 -b input.bam > mapped.bam
samtools sort mapped.bam > sorted.bam
```

Then, run plotCoverage:
```
plotCoverage -b sorted_sample1.bam  sorted_sample2.bam \
--plotFile plotCoverage.sample1ANDsample2_mapping_mouse.no_duplicates \
--plotTitle "Coverage after hybrid mapping" \
--ignoreDuplicates \
--numberOfProcessors 8 \
--labels sample1 sample2
```

### Some statistics with Covtobed
> A tool to generate BED coverage tracks from BAM files (![covtobed.sh](https://github.com/mvie-git/ISA-linked-reads/blob/main/scripts/MOLECULE-first/coverage/covtobed/covtobed.sh))

