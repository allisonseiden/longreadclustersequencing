# samtools adding read groups to files
samtools addreplacerg -r ID:ID number.1 -r LIB:hg38 -r SM:ID -o ID_edit.bam <bam_file>

# samtools making sure files are same size
samtools view file_name | wc -l

# make index files with python script, use subprocess? before calling whatshap
# for every individual, make their index file first
