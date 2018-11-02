#BSUB -W 3:00
#BSUB -q premium
#BUSB -n 2
#BSUB -R "rusage[mem=50000]"
#BSUB -P acc_schade01a
#BSUB -J whatshap_b3_11_01
#BSUB -m mothra
#BSUB -o whatshap_b3_11_01.stdout
#BSUB -e whatshap_b3_11_01.stderr


# submit with this command: 
# cd /hpc/users/richtf01/longreadclustersequencing/phasing
# for i in {1..200}; do echo $i; bsub < whatshap_bsub.sh; done

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch3/
python3 ~/longreadclustersequencing/phasing/illumina_whatshap_int1.py --batch 3
