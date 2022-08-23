# Extract list of barcodes (from mouse molecule file) associated with only one molecule
sed 1d 149_mapping_mouse.molecules.chimeric.bed | awk '{print $4}' | sort | uniq -c | sort -n | awk '$1==1 {print $2}' > 149_bx_one_molecule.txt

# Subset 
grep -f 149_bx_one_molecule.txt 149_mapping_mouse.molecules.chimeric.bed
grep -f 149_bx_one_molecule.txt 149_mapping_AAV.molecules.chimeric.bed

# Add prefix to barcodes
sed -i -e 's/^/BX:Z:/' 149_bx_one_molecule.txt

