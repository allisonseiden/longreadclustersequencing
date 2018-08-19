#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Split VCF by chromosome.

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
IlluminaWhatshapVCFs/Batch3/
python3 ~/longreadclustersequencing/phasing/illumina_whatshap_int1.py

deactivate
"""

import pandas as pd
import subprocess as sp
import multiprocessing as mp

""" Script 1/2 to run Whatshap with indels flag using Illumina data for first
    5 patient IDs, uses multiprocessing """

df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/' +
                   'Illumina_crams_batch3.txt', names=['ID'])
patientID = df['ID'].tolist()

# patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389"]
patientID = ['CG0000-1789']  # this is 1-00004


def illumina_whatshap(ID):
    """Run whatshap for illumina data."""
    print(ID)
    mkdir = 'mkdir -p ' + ID
    sp.call(mkdir, shell=True)
    cd = ('cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
          'IlluminaWhatshapVCFs/Batch1/' + ID)
    sp.call(cd, shell=True)
    for i in range(1, 23):
        # full VCF:
        # vcf_filename = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
        #                 'PacbioProject/GMKF_TrioVCFs/' + ID + '_trio.vcf.gz')
        # per chrom VCF:
        fam_id = '1-00004'
        vcf_filename = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                        'PacbioProject/GMKF_TrioVCFs/{}/Illumina_WGS_{}' +
                        '_chr{}.vcf.gz').format(fam_id, fam_id, i)
        bam_filename = ('/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_' +
                        '2017/CRAM/Batch3/' + ID + '.cram')
        command = ('time whatshap phase --sample=' + ID +
                   ' --ignore-read-groups' + ' --reference ' +
                   '/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa --indels ' +
                   '-o ' + ID + '/' + fam_id + '_chr' +
                   str(i) + '_phased.vcf ' + vcf_filename + ' ' + bam_filename)
        print(command)
        sp.call(command, shell=True)
        print('======Sucessfully ran whatshap for ' + ID + ' on chr ' + str(i))
        break
    sp.call('cd ..', shell=True)


if __name__ == '__main__':
    pool = mp.Pool(processes=5)
    pool.map(illumina_whatshap, patientID)
