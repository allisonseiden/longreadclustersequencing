#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Split VCF by sample.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-10
:Copyright: 2018, Allison Seiden
:License: CC BY-SA


module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/GMKF_TrioVCFs/
python3 ~/longreadclustersequencing/phasing/split_trios.py

Once done, check that last line in every file is the same

"""


import pandas as pd
import subprocess as sp
import multiprocessing as mp
import os


# batch_prefix = '/hpc/users/seidea02/longreadclustersequencing/data/Illumina_'
# batch1_df = pd.read_table(batch_prefix + 'crams_batch1.txt', names=['ID']);
# batch1 = batch1_df['ID'].tolist();
#
# batch2_df = pd.read_table(batch_prefix + 'crams_batch2.txt', names=['ID']);
# batch2 = batch2_df['ID'].tolist();
#
# batch3_df = pd.read_table(batch_prefix + 'crams_batch3.txt', names=['ID']);
# batch3 = batch3_df['ID'].tolist();


def split_vcf(num):
    """Split VCF into trios to avoid whatshap memory issues."""
    kiddo = trio_df.loc[num, 'Child'];
    dad = trio_df.loc[num, 'Father'];
    mom = trio_df.loc[num, 'Mother'];
    print(kiddo)
    vcf_exists = os.path.exists(kiddo + '_trio.vcf.gz')
    tbx_exists = os.path.exists(kiddo + '_trio.vcf.gz.tbi')
    # first check if VCF created but not indexed (likely due to crash)
    # DO NOT use if running on multiple screens
    # if vcf_exists and not tbx_exists:
    #     print('Deleting and restarting this file: ' + kiddo + '_trio.vcf.gz')
    #     os.remove(kiddo + '_trio.vcf.gz')
    # if BOTH VCF and tabix were created, skip and move on to next one
    if vcf_exists and tbx_exists:
        print('Already done with ' + kiddo)
        return kiddo
    command = 'bcftools view -s ' + kiddo + ',' + mom + ',' + dad;
    command += ' -O z -o ' + kiddo + '_trio.vcf.gz /sc/orga/projects/';
    command += 'chdiTrios/GMKF_WGS_Trios_Dec_2017/';
    command += 'GMKF_Seidman_CongenitalHeartDisease_WGS.vcf.gz';
    index = 'tabix -p vcf ' + kiddo + '_trio.vcf.gz';
    """
    cd = ('cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
          'PacbioProject/');
    if kiddo in batch1:
        sp.call(cd + 'GMKF_TrioVCFs/Batch1', shell=True);
    elif kiddo in batch2:
        sp.call(cd + 'GMKF_TrioVCFs/Batch2', shell=True);
    else:
        sp.call(cd + 'GMKF_TrioVCFs/Batch3', shell=True);
    sp.call(cd, shell=True);
    """
    sp.call(command, shell=True);
    sp.call(index, shell=True);
    return kiddo


if __name__ == '__main__':
    trio_df = pd.read_table('/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_' +
                            '2017/GMKF_Seidman_CongenitalHeartDisease_WGS.ped',
                            names=['Fam_ID', 'Child', 'Father', 'Mother'],
                            sep='\t', comment='#');
    length = trio_df.shape[0];
    pool = mp.Pool(processes=7);
    pool.map(split_vcf, range(length));
