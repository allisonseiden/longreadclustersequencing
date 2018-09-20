#BSUB -W 4:00
#BSUB -q alloc
#BUSB -n 2
#BSUB -R "rusage[mem=50000]"
#BSUB -P acc_chdiTrios
#BSUB -J whatshap_b1_t10
#BSUB -m mothra
#BSUB -o whatshap_b1_t10.stdout
#BSUB -e whatshap_b1_t10.stderr


# submit with this command:
# bsub < whatshap_bsub.sh

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch3/
python3 ~/longreadclustersequencing/phasing/illumina_whatshap_int1.py 3
