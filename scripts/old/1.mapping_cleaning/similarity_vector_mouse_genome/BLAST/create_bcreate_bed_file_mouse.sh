# Remove first line
sed '1d' blast_output_vector_mouse_vector_elements.txt > file.tmp

# Extract specific columns (mouse chr, start, end positions)
awk '{print $2 "\t" $9 "\t" $10}' file.tmp > blast_output_mouse_filter.tmp

# Remove double quotes in the chr column
cat blast_output_mouse_filter.tmp | tr -d '"' > blast_output_mouse_filter.bed

# Clean
rm *.tmp