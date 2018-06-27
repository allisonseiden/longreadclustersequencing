# how to remove non-ascii characters from bam
sed 's/[^\x00-\x7F]//g' ID_edit.bam > ID_edit.sam

# double checking there are no non-ascii characters
grep -P -n --color="auto" "[^\x00-\x7F]" ID_edit.sam
