**Chimeric-first** strategy consists in identifying chimeric reads then retrieve all reads sharing the same barcode to provide long-range information on the integration site. There is three steps:
1. Identify chimeric reads by using specific tools for LAM-PCR or TES integration site analysis
2. Retrieve all the barcodes associated with the chimeric reads identified to reconstruct the initial DNA molecule
3. Map all the reads gathered to provide more information around the chimeric read identified


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

## Integration site detection tools
> Tools have been developed to identify exact positions of human and virus breakpoints of integration events. Then it can annotate the integration events with the corresponding affected human genes. In our case, the goal is to study the integration events after injection of a vector in a mouse to treat a disease. No amplification and/or hybridization of ITRs/viral sequences yielding short reads have been applied.

After benchmarking existing tools, here is the list of selected tools:
- ![Virus-clip](https://github.com/mvie-git/ISA-linked-reads/tree/main/scripts/CHIMERIC-first/integration_site_detection_tools/virus_clip)
- ![GENE-IS]()
- 







