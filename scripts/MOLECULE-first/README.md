# ISA: Molecule-first strategy

This strategy consists in gathering all reads sharing the same barcode into molecules then identify chimeric molecules.
There is three steps:
1. Gather all reads sharing the same barcode to reconstruct the initial DNA molecule
2. Map all the molecules against the hybrid genome
3. Identify chimeric molecules meaning with reads mapped to both genomes


## Generate TELL-Seq linked-reads from sequencing data

> The Tell-Read module of the TELL-Seq software demultiplex and correct error barcoded reads from the initial raw sequencing data.

The Tell-Seq pipeline (`Tell-Read` module) takes as input the initial raw sequencing FASTQ data: 
- **demultiplexing**: they are further demultiplexed into sample separated I1 (index 1) reads, R1 reads, R2 reads, based on I2 (index 2) reads. I1 reads are the TELL_seq barcode sequences. For each sequencing library construction, a set of unique barcode sequences was randomly chosen from a 2.4 billion-barcode pool. These sample-demultiplexed FASTQ files are saved as the raw data output files. 
- **adaptor trimming**: adapter sequences are then trimmed using *cutadapt* utility.
- **error barcoded reads correction**: the unique barcodes associated with only one read are most likely caused by sequencing errors in the barcode. These barcodes are first identified if they are 1-base mismatches with one of the barcode associated with multiple reads, and then error-corrected. Barcodes with errors after this step are filtered out. The erroneous barcodes along with their associated reads are removed and excluded from the rest of analyses. The remainin R1 and R2 reads, along with their associated I1 reads (barcodes) are the TELL_seq linked reads. They are the input for downstream analysis

To run Tell-Read module:
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


## Molecule-first



