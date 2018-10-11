#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Run whatshap with indels per chromosome.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-19
:Copyright: 2018, Allison Seiden
:License: CC BY-SA


cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch2/
python3 ~/longreadclustersequencing/phasing/illumina_whatshap_int1.py --batch 2

CG0011-0730
# once done, need to confirm that last line is the same in every file
# e.g.,  zcat GMKF_TrioVCFs/1-00004/Illumina_WGS_1-00004_chr4.vcf.gz | tail
# and IlluminaWhatshapVCFs/Batch3/CG0000-1789/1-00004_chr4_phased.vcf

deactivate
"""

import os
import subprocess as sp
# import multiprocessing as mp
import argparse
# from functools import partial

from utils import get_trio_df, get_batch_pt_ids


""" Script 1/2 to run Whatshap with indels flag using Illumina data for first
    5 patient IDs, uses multiprocessing """


def illumina_whatshap_per_chrom(ID, batch_ct):
    """Run whatshap for illumina data by CHROMOSOME."""
    # obtain the family ID
    trio_df = get_trio_df()
    fam_id = trio_df.loc[trio_df.Child == ID]['Fam_ID'].to_string(index=False)
    print('Starting {} (family {})'.format(ID, fam_id))
    # make and enter directories named after the DNA ID
    mkdir = 'mkdir -p ' + ID
    sp.call(mkdir, shell=True)
    cd = ('cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
          'IlluminaWhatshapVCFs/Batch{}/{}').format(batch_ct, ID)
    print(cd)
    sp.call(cd, shell=True)
    for i in range(1, 23):
        vcf_filename = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                        'PacbioProject/GMKF_TrioVCFs/{}/Illumina_WGS_{}' +
                        '_chr{}.vcf.gz').format(fam_id, fam_id, i)
        if not os.path.exists(vcf_filename):
            print('chr{} from {} not ready yet'.format(i, ID))
            continue
        bam_filename = ('/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_' +
                        '2017/CRAM/Batch{}/{}.cram').format(batch_ct, ID)
        command = ('time whatshap phase --sample=' + ID +
                   ' --ignore-read-groups --reference ' +
                   '/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa --indels ' +
                   '-o ' + ID + '/' +
                   fam_id + '_chr' +
                   str(i) + '_phased.vcf ' + vcf_filename + ' ' + bam_filename)
        if os.path.exists('{}/{}_chr{}_phased.vcf'.format(ID, fam_id, i)):
            print('chr{} from {} already run'.format(i, ID))
            # check file size (only for cleaning)
            out_f = '{}/{}_chr{}_phased.vcf'.format(ID, fam_id, i)
            print(os.stat(out_f).st_size)
            # delete if less than 100 Mb, suggests it did not complete
            # if os.stat(out_f).st_size == 0:
            #     print('deleting unfinishing phasing for ' + out_f)
            #     # os.remove(out_f)
            # else:
            continue
        print(command)
        sp.call(command, shell=True)
        print('======Sucessfully ran whatshap for ' + ID + ' on chr ' + str(i))
    sp.call('cd ..', shell=True)


if __name__ == '__main__':
    # pool = mp.Pool(processes=2)
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch', default=1, type=int,
                        choices=[1, 2, 3], help='Pick the batch')
    args = parser.parse_args()
    batch_ct = args.batch
    cd = ('cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
          'IlluminaWhatshapVCFs/Batch' + str(batch_ct))
    sp.call(cd, shell=True)
    patientID_list = get_batch_pt_ids(batch_ct)  # [:10]
    # whatshap_partial = partial(illumina_whatshap_per_chrom,batch_ct=batch_ct)
    # done_ids = pool.map(whatshap_partial, patientID)
    done_ids = []
    for patientID in patientID_list:
        illumina_whatshap_per_chrom(patientID, batch_ct)
    print(done_ids)
