# notes for Friday June 8th

# running whatshap
whatshap phase --ped test.ped --indels --reference
/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa --distrust-genotypes --include-homozygous -o
phasedFile.vcf ~/www/PacbioProject/DNV_calls/VCF/ID ~/www/PacbioProject/PacbioHg38Bams/ID
~/www/PacbioProject/IlluminaHg38Bams/ID-mom ~/www/PacbioProject/IlluminaHg38Bams/ID-dad

# using picard
module load picard
java -jar $PICARD AddOrReplaceReadGroups INPUT=~/www/PacbioProject/PacbioHg38Bams/ID
OUTPUT=~/longreadclustersequencing/test_bam.bam RGLB=lib1 RGPL=pacbio RGPU=unit1 RGSM=id number
VALIDATION_STRINGENCY=SILENT (possibly what will make picard ignore read quality)

java -jar $PICARD BuildBamIndex INPUT=bam_file


# samtools
samtools view file__name | wc -l
