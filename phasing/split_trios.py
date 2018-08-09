import pandas as pd
import subprocess as sp
import multiprocessing as mp


batch1_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/Illumina_crams_batch1.txt', names=['ID']);
batch1 = batch1_df['ID'].tolist();

batch2_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/Illumina_crams_batch2.txt', names=['ID']);
batch2 = batch2_df['ID'].tolist();

batch3_df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/Illumina_crams_batch3.txt', names=['ID']);
batch3 = batch3_df['ID'].tolist();

trio_df = pd.read_table('/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_2017/GMKF_Seidman_CongenitalHeartDisease_WGS.ped',
                        names = ['Fam_ID', 'Child', 'Father', 'Mother'], sep='\t', comment='#');

length = trio_df.shape[0];
def split_vcf(num):
    kiddo = trio_df.loc[num, 'Child'];
    dad = trio_df.loc[num, 'Father'];
    mom = trio_df.loc[num, 'Mother'];
    command = 'bcftools view -s ' + kiddo + ',' + mom + ',' + dad;
    command += '-O z -o ' + kiddo + '_trio.vcf.gz /sc/orga/projects/';
    command += 'chdiTrios/GMKF_WGS_Trios_Dec_2017/';
    command += 'GMKF_Seidman_CongenitalHeartDisease_WGS.vcf.gz';
    index = 'tabix -p vcf ' + kiddo + '_trio.vcf.gz';
    if kiddo in batch1:
        cd = 'cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/';
        cd += 'GMKF_TrioVCFs/Batch1';
    elif kiddo in batch2:
        cd = 'cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/';
        cd += 'GMKF_TrioVCFs/Batch2';
    else:
        cd = 'cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/';
        cd += 'GMKF_TrioVCFs/Batch3';
    sp.call(cd, shell=True);
    sp.call(command, shell=True);
    sp.call(index, shell=True);


if __name__ == '__main__':
    pool = mp.Pool(processes=5);
    pool.map(split_vcf, range(length));
