workdir="/home/mallaury/Bureau/phd/project/GSDIa_louisa/ISA"

while read sample; do

	echo "*** Start $sample ***"
	
	input="$workdir/4.identify_chimeric_molecules/output/${sample}/${sample}_mapping_mouse.molecules.chimeric.bed"

	# Histogram on number of molecules per barcode
	sed 1d $input | awk '{print $4}' | sort | uniq -c | awk '{print $1}' | awk '{print '${sample}' "\t" $0}' > $workdir/visualization/r_ggplot_boxplot/input/${sample}.csv


done < sample_list.txt


# Gather all data
cat $workdir/visualization/r_ggplot_boxplot/input/* > $workdir/visualization/r_ggplot_boxplot/input/all.csv