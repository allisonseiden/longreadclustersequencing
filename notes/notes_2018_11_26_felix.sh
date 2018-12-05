

# Problematic IDs for Batch 1 (all had KeyErrors)
['1-06459', '1-11584', '1-14377', '1-10377', '1-08634', '1-00101', '1-03142', '1-13125', '1-14851', '1-11319', '1-15514', '1-00030', '1-13588', '1-15150']
# incomplete dataframes
# 1-00186_dataframe_incomplete.txt, 1-00599_dataframe_incomplete.txt, 1-08416_dataframe_incomplete.txt, 1-15267_dataframe_incomplete.txt, 1-15358_dataframe_incomplete.txt


diff <(wc -l gmkf2/*) <(wc -l gmkf2_old/* | sed 's/_old//g')

# TESTING sorting hat

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix/0.2.6 bedtools/2.27.1
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate

# only run once on login nodes:
# pip install sorting_hat
sorting_hat

cd /hpc/users/richtf01/longreadclustersequencing/data/gmkf2

FA='/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa'
RMSK='/hpc/users/richtf01/longreadclustersequencing/data/repeat_masker.txt'
time sorting_hat --bed 1-00004_dnv.bed \
                 --fasta $FA \
                 --repeat $RMSK \
                 --output 'sorting_hat_out/1-00004_dnv.txt'


