module purge
module load python/3.5.0 py_packages/3.5 tabix samtools/1.8

virtualenv venv_phasing
# confirm you see this as the first line:
# Using base prefix '/hpc/packages/minerva-common/python/3.5.0'
source venv_phasing/bin/activate
pip install --upgrade pip
pip install whatshap


# check output file sizes and completion with whatshap_output_check.py
# if per-trio VCF was split, move to done_2018_09_26 folder
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/GMKF_TrioVCFs/
# done_2018_09_26 done_2018_11_03
time tar -zcvf done_2018_11_03_archive.tar.gz done_2018_11_03

# if phasing is done, move split chrom to folder so it can be archived
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/GMKF_TrioVCFs/
# split_chr_done_2018_09_26
time tar -zcvf split_chr_done_2018_11_03_archive.tar.gz split_chr_done_2018_11_03

# must be on one of the internal login nodes: login1, login2, minerva4
# if not on these nodes then ssh in from minerva2
ssh login1
time dsmc archive -se=richtf01 done_2018_11_03_archive.tar.gz
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/GMKF_TrioVCFs/
for file in $(find split_chr_done_2018_11_03 -name '*vcf*'); do
    echo $file
    time dsmc archive -se=richtf01 $file
    rm $file
done

dsmc q archive -se=richtf01 /sc/orga/ -sub=yes | grep "split_chr_done_2018_09_26"

# Once confident that phasing is complete for current dataset,
# remove line if 0|1 1|0 0/1 or 1/1 are NOT present (i.e., only 0/0 or ./.)

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/IlluminaWhatshapVCFs/Batch2/CG0011-3774/
time grep "^#\|/1\||1\|1|" 1-04266_chr14_phased.vcf
time grep "|1\|1|" 1-04266_chr14_phased.vcf
time grep "|1\|1|" 1-04266_chr22_phased.vcf | wc -l


time whatshap phase --sample=CG0007-7707 --ignore-read-groups --reference /sc/orga/projects/chdiTrios/Felix/db
s/hg38.fa --indels -o CG0007-7707/1-03560_chr17_phased.vcf /sc/orga/projects/chdiTrios/WGS_Combined_2017/Pacbi
oProject/GMKF_TrioVCFs/1-03560/Illumina_WGS_1-03560_chr17.vcf.gz /sc/orga/projects/chdiTrios/GMKF_WGS_Trios_De
c_2017/CRAM/Batch3/CG0007-7707.cram
