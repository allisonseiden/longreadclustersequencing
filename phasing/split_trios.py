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

Once done, check that last line in every file is the same.

There are 16 that still need to be split, but because of Minerva space
issues I am moving the 400 already split to a separate directory
(GMKF_TrioVCFs/done_2018_09_26) that will be archived

"""


import pandas as pd
import subprocess as sp
import multiprocessing as mp
import os


# batch_prefix = '/hpc/users/seidea02/longreadclustersequencing/data/Illumina_'
# batch1_df = pd.read_table(batch_prefix + 'crams_batch1.txt', names=['ID'])
# batch1 = batch1_df['ID'].tolist()
#
# batch2_df = pd.read_table(batch_prefix + 'crams_batch2.txt', names=['ID'])
# batch2 = batch2_df['ID'].tolist()
#
# batch3_df = pd.read_table(batch_prefix + 'crams_batch3.txt', names=['ID'])
# batch3 = batch3_df['ID'].tolist()


final_set = ['CG0009-7084', 'CG0015-5533', 'CG0017-4546', 'CG0020-9963',
             'CG0021-6581', 'CG0020-1859', 'CG0020-8872', 'CG0023-5663',
             'CG0021-2288', 'CG0023-0457', 'CG0021-3351', 'CG0018-5622',
             'CG0022-6769', 'CG0019-8485', 'CG0021-6490', 'CG0020-0171']


def split_vcf(num):
    """Split VCF into trios to avoid whatshap memory issues."""
    kiddo = trio_df.loc[num, 'Child']
    dad = trio_df.loc[num, 'Father']
    mom = trio_df.loc[num, 'Mother']
    vcf_exists = os.path.exists(kiddo + '_trio.vcf.gz')
    vcf_done = os.path.exists('done_2018_09_26/' + kiddo + '_trio.vcf.gz')
    """
    tbx_exists = os.path.exists(kiddo + '_trio.vcf.gz.tbi')
    # check if VCF created but not indexed (e.g., due to crash)
    # DO NOT use if running on multiple screens
    if vcf_exists and not tbx_exists:
        print('Deleting and restarting this file: ' + kiddo + '_trio.vcf.gz')
        os.remove(kiddo + '_trio.vcf.gz')
    # """
    # if VCF was created, skip and move on to next one
    if vcf_exists or vcf_done:
        # print('Already done with ' + kiddo)
        return 'not_rerun_' + kiddo
    if kiddo not in final_set:
        # need this if since since vcf_done and vcf_exists were archived
        return 'not_rerun_' + kiddo
    print('Starting ' + kiddo)
    command = 'bcftools view -s ' + kiddo + ',' + mom + ',' + dad
    command += ' -O z -o ' + kiddo + '_trio.vcf.gz /sc/orga/projects/'
    command += 'chdiTrios/GMKF_WGS_Trios_Dec_2017/'
    command += 'GMKF_Seidman_CongenitalHeartDisease_WGS.vcf.gz'
    index = 'tabix -p vcf ' + kiddo + '_trio.vcf.gz'
    """
    cd = ('cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
          'PacbioProject/')
    if kiddo in batch1:
        sp.call(cd + 'GMKF_TrioVCFs/Batch1', shell=True)
    elif kiddo in batch2:
        sp.call(cd + 'GMKF_TrioVCFs/Batch2', shell=True)
    else:
        sp.call(cd + 'GMKF_TrioVCFs/Batch3', shell=True)
    sp.call(cd, shell=True)
    """
    print(command)
    sp.call(command, shell=True)
    print(index)
    sp.call(index, shell=True)
    return kiddo


if __name__ == '__main__':
    trio_df = pd.read_table('/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_' +
                            '2017/GMKF_Seidman_CongenitalHeartDisease_WGS.ped',
                            names=['Fam_ID', 'Child', 'Father', 'Mother'],
                            sep='\t', comment='#')
    length = trio_df.shape[0]
    pool = mp.Pool(processes=3)
    kiddo_done = pool.map(split_vcf, range(length))
    print('Total number of probands: {}'.format(len(kiddo_done)))
    kiddo_previously_done = [i for i in kiddo_done if 'not_rerun_' in i]
    print('Probands completed: {}'.format(len(kiddo_previously_done)))
    # print('Probands to do:')
    # print([i for i in kiddo_done if 'not_rerun_' not in i])
