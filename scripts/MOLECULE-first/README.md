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

## BWA: Map TELL-seq linked reads against the hybrid genome
> BWA is a software package for mapping low-divergent sequences against a large reference genome. BWA-MEM is better for longer sequences ranged from 70bp to 1Mbp.

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




