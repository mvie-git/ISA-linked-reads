while read sample; do

	echo "*** Start $sample ***"

	# Set up parameters
	output_virus_clip="result/${sample}/virus_clip.out"
	BED_mouse_similar_regions="input/mouse_blacklist_regions_similar.bed"
	BED_vector_similar_regions="input/vector_blacklist_regions_similar.sorted.bed"

	# Create the output directory
	mkdir -p result/$sample/remove_FP

	# # DONE! But before the BED file need to be sorted by chromosome and then by start position
	# sed 's/chr//g' mouse_blacklist_regions_similar.original.txt | sed 1d > mouse_blacklist_regions_similar.bed

	# If the left-element of the soft-clipped read is a mouse, extract the integration loci and calculate the length of the clipped read to guess the end of the mouse integration loci (part 1)
	echo "*** Starting the extraction of the integration loci if the left-element of the soft-clipped read is a mouse ***"
	awk '$1=="mouse" {print $2 "\t" $3 "\t" $4}' $output_virus_clip | awk '{print $0"\t"length($3)}' | awk '$5=$2+$4' | awk '{print $1 "\t" $2 "\t" $5}' > mouse_integration_loci_left.bed
	echo "*** Extract the left-element of soft-clipped reads done ***"

	# If the right-element of the soft-clipped read is a mouse, extract the integration loci and calculate the length of the clipped read to guess the end of the mouse integration loci (part 2)
	echo "*** Starting the extraction of the integration loci if the right-element of the soft-clipped read is a mouse ***"
	awk '$5=="mouse" {print $6 "\t" $7 "\t" $8}' $output_virus_clip | awk '{print $0"\t"length($3)}' | awk '$5=$2+$4' | awk '{print $1 "\t" $2 "\t" $5}' > mouse_integration_loci_right.bed
	echo "*** Extract the right-element of soft-clipped reads done ***"

	# Concatenate the left and right clipped-reads
	echo "*** Starting the concatenation of the left and right clipped-reads ***"
	cat mouse_integration_loci_left.bed mouse_integration_loci_right.bed > mouse_integration_loci.bed
	echo "*** Concatenate the soft-clipped done ***"

	# Sort BED file which contained the clipped-reads by chromosome and then by start position
	echo "*** Starting the sorting of the BED file of clipped-reads by start position ***"
	sort -k1,1 -k2,2n mouse_integration_loci.bed > result/$sample/remove_FP/mouse_integration_loci.sorted.bed
	echo "*** Sort BED file done ***"

	# Bedtools intersect: check if any of the features in the two genomic sets “overlap” with one another

	# Attention! Before executing the bedtools intersect command, be sure that the file is TAB delimited
	echo "*** Format the BED file to be TAB-delimited for bedtools command ***"
	awk '{print $1 "\t" $2 "\t" $3}' result/$sample/remove_FP/mouse_integration_loci.sorted.bed > result/$sample/remove_FP/mouse_integration_loci.sorted.tab.bed
	echo "*** TAB-delimited BED file done"

	## Options:
	### -a: Each feature in A is compared to B in search of overlaps (BAM/BED/GFF/VCF file)
	### -b: One or more BAM/BED/GFF/VCF file
	### -v: Only report those entries in A that have no overlap in B (here only report integration loci that are not in similar regions between the two references)
	echo "*** Starting the execution of the bedtools intersect tool ***"
	bedtools intersect -v -a result/$sample/remove_FP/mouse_integration_loci.sorted.tab.bed -b $BED_mouse_similar_regions > result/$sample/remove_FP/bedtools_intersect_mouse.bed
	echo "*** Bedtools intersect is done ***"

	# Comm: compare two sorted files line by line and write to standard output; the lines that are common and the lines that are unique

	# Attention! Before executing the comm command, be sure that the two files you compare are sorted
	# Sort BED file which contained the clipped-reads by chromosome and then by start position
	echo "*** Starting the sorting of the BED file of integration loci by chromosome then by the start position ***"
	sort -k1,1 -k2,2n result/$sample/remove_FP/mouse_integration_loci.sorted.tab.bed > result/$sample/remove_FP/mouse_integration_loci.sorted.tab.sorted.bed
	echo "*** Sort BED files done ***"

	## Options:
	### -2: suppress second column (lines unique to second file)
	### -3: suppress third column (lines unique to both files)
	# Keep only the first column which contained only the integration loci in the blacklist region
	echo "*** Extract the integration loci which are in a blacklist region ***"
	comm -23 result/$sample/remove_FP/mouse_integration_loci.sorted.tab.sorted.bed result/$sample/remove_FP/bedtools_intersect_mouse.bed > result/$sample/remove_FP/bedtools_intersect_mouse_FP.bed
	echo "*** BED file contained integration loci of false positives created ***"

	# Format the integration loci in a text file to allow the grep command to remove false positives from the soft-clipped reads file
	echo "*** Format the integration loci file to use grep command then ***"
	awk '{print $1 "\t" $2}' result/$sample/remove_FP/bedtools_intersect_mouse_FP.bed > result/$sample/remove_FP/mouse_FP_to_remove.txt
	echo "*** Created text file with false positives done ***"

	# Filter out the false positives from the output of Virus-clip software
	# Output: 1) Virus clip output file without false positives 2) Virus clip output contained only false positives
	echo "*** Starting to filter out false positives from the output of Virus Clip software ***"
	sed 1d $output_virus_clip | grep -f result/$sample/remove_FP/mouse_FP_to_remove.txt > result/$sample/remove_FP/virus_clip.FP_mouse.out
	sed 1d $output_virus_clip | grep -v -f result/$sample/remove_FP/mouse_FP_to_remove.txt > result/$sample/remove_FP/virus_clip.filter_mouse.out
	echo "*** Filtering out false positives is done ***"

	# Clean intermediate files
	echo "*** Cleaning intermediate files ***"
	rm mouse_integration_loci_left.bed
	rm mouse_integration_loci_right.bed
	rm mouse_integration_loci.bed
	rm result/$sample/remove_FP/mouse_integration_loci.sorted.bed
	rm result/$sample/remove_FP/mouse_integration_loci.sorted.tab.bed
	rm result/$sample/remove_FP/mouse_FP_to_remove.txt


	#------------------------------------------------------------------------------------------------------------------------------------#

	# Do the same thing with soft-elements mapped to the vector

# # DONE! But before the BED file need to be sorted by chromosome and then by start position
	# sed 1d vector_blacklist_regions_similar.txt | sort -k1,1 -k2,2n > vector_blacklist_regions_similar.sorted.bed

	# If the left-element of the soft-clipped read is a vector, extract the integration loci and calculate the length of the clipped read to guess the end of the vector integration loci (part 1)
	echo "*** Starting the extraction of the integration loci if the left-element of the soft-clipped read is a vector ***"
	awk '$1=="Virus" {print $2 "\t" $3 "\t" $4}' $output_virus_clip | awk '{print $0"\t"length($3)}' | awk '$5=$2+$4' | awk '{print $1 "\t" $2 "\t" $5}' | awk '$1="1"' > vector_integration_loci_left.bed
	echo "*** Extract the left-element of soft-clipped reads done ***"

	# If the right-element of the soft-clipped read is a vector, extract the integration loci and calculate the length of the clipped read to guess the end of the vector integration loci (part 2)
	echo "*** Starting the extraction of the integration loci if the right-element of the soft-clipped read is a vector ***"
	awk '$5=="Virus" {print $6 "\t" $7 "\t" $8}' $output_virus_clip | awk '{print $0"\t"length($3)}' | awk '$5=$2+$4' | awk '{print $1 "\t" $2 "\t" $5}' | awk '$1="1"' > vector_integration_loci_right.bed
	echo "*** Extract the right-element of soft-clipped reads done ***"

	# Concatenate the left and right clipped-reads
	echo "*** Starting the concatenation of the left and right clipped-reads ***"
	cat vector_integration_loci_left.bed vector_integration_loci_right.bed > vector_integration_loci.bed
	echo "*** Concatenate the soft-clipped done ***"

	# Sort BED file which contained the clipped-reads by chromosome and then by start position
	echo "*** Starting the sorting of the BED file of clipped-reads by start position ***"
	sort -k1,1 -k2,2n vector_integration_loci.bed > result/$sample/remove_FP/vector_integration_loci.sorted.bed
	echo "*** Sort BED file done ***"

	# Bedtools intersect: check if any of the features in the two genomic sets “overlap” with one another

	# Attention! Before executing the bedtools intersect command, be sure that the file is TAB delimited
	echo "*** Format the BED file to be TAB-delimited for bedtools command ***"
	awk '{print $1 "\t" $2 "\t" $3}' result/$sample/remove_FP/vector_integration_loci.sorted.bed > result/$sample/remove_FP/vector_integration_loci.sorted.tab.bed
	echo "*** TAB-delimited BED file done"

	## Options:
	### -a: Each feature in A is compared to B in search of overlaps (BAM/BED/GFF/VCF file)
	### -b: One or more BAM/BED/GFF/VCF file
	### -v: Only report those entries in A that have no overlap in B (here only report integration loci that are not in similar regions between the two references)
	echo "*** Starting the execution of the bedtools intersect tool ***"
	bedtools intersect -v -a result/$sample/remove_FP/vector_integration_loci.sorted.tab.bed -b $BED_vector_similar_regions > result/$sample/remove_FP/bedtools_intersect_vector.bed
	echo "*** Bedtools intersect is done ***"

	# Comm: compare two sorted files line by line and write to standard output; the lines that are common and the lines that are unique
	## Options:
	### -2: suppress second column (lines unique to second file)
	### -3: suppress third column (lines unique to both files)
	# Keep only the first column which contained only the integration loci in the blacklist region
	echo "*** Extract the integration loci which are in a blacklist region ***"
	comm -23 result/$sample/remove_FP/vector_integration_loci.sorted.tab.bed result/$sample/remove_FP/bedtools_intersect_vector.bed > result/$sample/remove_FP/bedtools_intersect_vector_FP.bed
	echo "*** BED file contained integration loci of false positives created ***"

	# Format the integration loci in a text file to allow the grep command to remove false positives from the soft-clipped reads file
	echo "*** Format the integration loci file to use grep command then ***"
	awk '$1="chrVirus"' result/$sample/remove_FP/bedtools_intersect_vector_FP.bed | awk '{print $1 "\t" $2}' > result/$sample/remove_FP/vector_FP_to_remove.txt
	echo "*** Created text file with false positives done ***"

	# Filter out the false positives from the output of Virus-clip software
	# Output: 1) Virus clip output file without false positives 2) Virus clip output contained only false positives
	echo "*** Starting to filter out false positives from the output of Virus Clip software ***"
	grep -f result/$sample/remove_FP/vector_FP_to_remove.txt result/$sample/remove_FP/virus_clip.filter_mouse.out > result/$sample/remove_FP/virus_clip.FP_vector.out
	grep -v -f result/$sample/remove_FP/vector_FP_to_remove.txt result/$sample/remove_FP/virus_clip.filter_mouse.out > result/$sample/remove_FP/virus_clip.filter_FP.out
	echo "*** Filtering out false positives is done ***"

	# Concatenate FP from mouse and vector
	cat result/$sample/remove_FP/virus_clip.FP_mouse.out result/$sample/remove_FP/virus_clip.FP_vector.out > result/$sample/remove_FP/virus_clip.FP.out

	# Clean intermediate files
	echo "*** Cleaning intermediate files ***"
	rm vector_integration_loci_left.bed
	rm vector_integration_loci_right.bed
	rm vector_integration_loci.bed
	rm result/$sample/remove_FP/vector_integration_loci.sorted.bed
	rm result/$sample/remove_FP/virus_clip.FP_vector.out
	rm result/$sample/remove_FP/virus_clip.FP_mouse.out
	rm result/$sample/remove_FP/vector_FP_to_remove.txt
	rm result/$sample/remove_FP/virus_clip.filter_mouse.out


	echo "*** Script is done. ***"

done < sample_list_remove_FP.txt