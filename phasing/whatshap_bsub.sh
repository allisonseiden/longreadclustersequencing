#BSUB -W 3:00
#BSUB -q premium
#BUSB -n 2
#BSUB -R "rusage[mem=50000]"
#BSUB -P acc_schade01a
#BSUB -J whatshap_b1_11_09
#BSUB -m mothra
#BSUB -o whatshap_b1_11_09.stdout
#BSUB -e whatshap_b1_11_09.stderr


# submit with this command: 
# cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
# IlluminaWhatshapVCFs
# cp ~/longreadclustersequencing/phasing/whatshap_bsub.sh .
# for i in {1..100}; do echo $i; bsub < whatshap_bsub.sh; done

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch3/
python3 ~/longreadclustersequencing/phasing/illumina_whatshap_int1.py --batch 3
