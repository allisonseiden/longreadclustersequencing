# commands for splitting by chromosome

module load tabix
cd into folder for current patient ID

# splitting
tabix -h ../ID.hg38.trio.vcf.gz chr_num > ID_chr_num.vcf

# compressing
bgzip ID_chr_num.vcf

# indexing
tabix -p vcf ID_chr_num.vcf.gz
