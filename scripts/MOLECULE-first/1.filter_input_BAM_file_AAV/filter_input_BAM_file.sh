# (Optional) Clean BAM file to filter out reads due to technical biases (unmapped, duplicates, multi-mapped) & extract AAV mapped reads

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

nb_threads=12

sambamba=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/bin/sambamba


## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Output directory
	mkdir -p $workdir/1.filter_input_BAM_file_AAV/${sample}

	# Input files
	input_BAM=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/output/${sample}/${sample}_mapping_mouseANDvector.sorted.chrVirus.bam

	# Filter out alignments which correspond to the FLAG option
	## Tool: samtools view
	## Options:
	### -b: output in the BAM format
	### -F: filter out alignments with all bits in INT present in the FLAG field
	### -@: number of threads
	
	# ### grep: select reads of interest
	# samtools view -@ $nb_threads $input_BAM | awk '$3=="chrVirus" {print $0}' > ${sample}_body.sam
	# samtools view -H $input_BAM > ${sample}_header.sam
	# cat ${sample}_header.sam ${sample}_body.sam > ${sample}_all.sam
	# samtools view -Sb ${sample}_all.bam
	# samtools sort ${sample}_all.bam -o $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.sortedchrVirus.mapped.bammapped.bam
	# samtools index $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.sortedchrVirus.mapped.bammapped.bam

	samtools view -b -F 4 $input_BAM > $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.bam

	# Identify and remove duplicate reads (PCR and optical duplicates) from the AAV mapping BAM file
	## Tool: sambamba markdup
	## Options:
	### -r: remove duplicates instead of just marking them
	### -t: number of threads
	### -p: show progressbar in the terminal
	$sambamba markdup -r -t $nb_threads -p $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.bam $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.bam

	### FLAG: 256 -> remove secondary alignments from the AAV mapping BAM file
	samtools view -@ $nb_threads -b -F 256 $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.bam > $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.no_secondary.bam 

	### FLAG: 2048 -> remove supplementary alignments from the AAV mapping BAM file
	samtools view -@ $nb_threads -b -F 2048 $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.no_secondary.bam > $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.no_secondary.no_supplementary.bam

	### Options: -q 1 -> remove reads with a mapping quality equal to zero (lower than 1)
	samtools view -b -q 1 $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.no_secondary.no_supplementary.bam > $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam


	# Sort and index the output filtered BAM file
	samtools sort -@ $nb_threads $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.bam -o $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.sorted.bam
	samtools index $workdir/1.filter_input_BAM_file_AAV/${sample}/${sample}_mapping_mouseANDvector.chrVirus.mapped.no_duplicates.sorted.bam


done < sample_list_1.txt
