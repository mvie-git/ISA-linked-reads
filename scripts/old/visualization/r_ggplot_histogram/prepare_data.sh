workdir="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA"

while read sample; do

	echo "*** Start $sample ***"
	
	input="$workdir/4.identify_chimeric_molecules/output/${sample}/${sample}_mapping_AAV.molecules.chimeric.bed"

	# Histogram on number of reads per molecule data
	sed 1d $input | awk '{print $5}' | awk '{print '${sample}' "\t" $0}' > $workdir/visualization/r_ggplot_histogram/input/${sample}.csv

	# # Histogram on length of molecules
	# sed 1d $input | awk '{print $6}' | awk '{print '${sample}' "\t" $0}' > $workdir/visualization/r_ggplot_histogram/input/${sample}.csv

done < sample_list.txt


# Gather all data
cat $workdir/visualization/r_ggplot_histogram/input/* > $workdir/visualization/r_ggplot_histogram/input/all.csv