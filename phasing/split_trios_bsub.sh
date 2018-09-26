#BSUB -W 4:00
#BSUB -q alloc
#BUSB -n 12
#BSUB -R "rusage[mem=50000]"
#BSUB -P acc_chdiTrios
#BSUB -J split_trios
#BSUB -m mothra
#BSUB -o split_trios.stdout
#BSUB -e split_trios.stderr


#  bsub < split_trios_bsub.sh 

module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/GMKF_TrioVCFs/
python3 ~/longreadclustersequencing/phasing/split_trios.py
