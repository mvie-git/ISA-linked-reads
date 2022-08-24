# (Optional) Clean BAM file to filter out reads due to technical biases (unmapped, duplicates, multi-mapped) & extract mouse mapped reads

workdir="/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/add_barcode_to_BAM"

nb_threads=12

mkdir -p $workdir/input/BAM

sambamba=/media/mallaury/MALLAURY1/PROJECT/GSDIa_ISA/hybrid_genome_mapping/bin/sambamba


## FOR EACH SAMPLE
while read sample; do

	echo "*** Start $sample ***"

	# Create an output directory
	mkdir -p ${sample}

	# Input files
	input_BAM=${workdir}/5.subset_BAM_file_with_readID/${sample}/${sample}_subset_BAM.sorted.bam

	# Filter out alignments which correspond to the FLAG option
	## Tool: samtools view
	## Options:
	### -b: output in the BAM format
	### -F: filter out alignments with all bits in INT present in the FLAG field
	### -@: number of threads
	
	### FLAG: 4 -> remove unmapped reads from the AAV mapping BAM file
	samtools view -@ $nb_threads -b -F 4 $input_BAM > ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.bam

	### grep: select reads of interest
	samtools view -@ $nb_threads ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.bam | awk '$3!="chrVirus" {print $0}' > ${sample}_body.sam
	samtools view -H $input_BAM > ${sample}_header.sam
	cat ${sample}_header.sam ${sample}_body.sam > ${sample}_all.sam
	samtools view -Sb ${sample}_all.sam > ${sample}_all.bam
	samtools sort ${sample}_all.bam -o ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.bam
	samtools index ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.bam

	# Clean
	rm ${sample}_body.sam ${sample}_header.sam ${sample}_all.sam ${sample}_all.bam

	# Identify and remove duplicate reads (PCR and optical duplicates) from the AAV mapping BAM file
	## Tool: sambamba markdup
	## Options:
	### -r: remove duplicates instead of just marking them
	### -t: number of threads
	### -p: show progressbar in the terminal
	$sambamba markdup -r -t $nb_threads -p ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.bam ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.bam

	# ### FLAG: 256 -> remove secondary alignments from the AAV mapping BAM file
	# samtools view -@ $nb_threads -b -F 256 ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.bam > ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.bam 

	# ### FLAG: 2048 -> remove supplementary alignments from the AAV mapping BAM file
	# samtools view -@ $nb_threads -b -F 2048 ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.bam > ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.bam

	# ### Options: -q 1 -> remove reads with a mapping quality equal to zero (lower than 1)
	# samtools view -b -q 1 ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.bam > ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam


	# # Sort and index the output filtered BAM file
	# samtools sort -@ $nb_threads ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.MAPQsup0.bam -o ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.MAPQsup0.sorted.bam
	# samtools index ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.no_secondary.no_supplementary.MAPQsup0.sorted.bam

	# Sort and index the output filtered BAM file
	samtools sort -@ $nb_threads ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.bam -o ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.sorted.bam
	samtools index ${workdir}/6.filter_input_BAM_file_mouse/${sample}/${sample}_subset_BAM.mapped.mouse.no_duplicates.sorted.bam

done < sample_list_6.txt
