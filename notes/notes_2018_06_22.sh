# running whatshap without pedigree mode

whatshap phase --sample=ID --ignore-read-groups
--reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa
--indels -o ID_chr22_phased.vcf
/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/DNV_calls/VCF/TrioVCF/ID/ID_chr22.vcf.gz
/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/PacbioHg38Bams/ID_edit_sorted.bam
