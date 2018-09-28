#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Split VCF by chromosome.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-19
:Copyright: 2018, Allison Seiden
:License: CC BY-SA


module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/GMKF_TrioVCFs/
python3 ~/longreadclustersequencing/phasing/chromosome_split.py

"""

import subprocess as sp
import os
# import multiprocessing as mp
from functools import partial
import gc

import pysam

from utils import get_trio_df

# Script to split VCF files into 22 separate VCF files by chromosome number

# patientID = ["1-01019", "1-03897", "1-04190", "1-04389", "1-04460",
#              "1-04537", "1-05443", "1-05673", "1-05846"];

# for ID in patientID:
#     cd = ("cd /hpc/users/seidea02/www/PacbioProject/DNV_calls/VCF/" +
#           "TrioVCF/" + ID);
#     subprocess.call(cd, shell=True);
#     for i in range(1, 23):
#         filename = ID + "_chr" + str(i) + ".vcf";
#         split = ("tabix -h /hpc/users/seidea02/www/PacbioProject/" +
#                  "DNV_calls/VCF/TrioVCF/" + ID + ".hg38.trio.vcf.gz chr" +
#                  str(i) + " > " + filename);
#         subprocess.call(split, shell=True);
#         compress = "bgzip " + filename;
#         subprocess.call(compress, shell=True);
#         index = "tabix -p vcf " + filename + ".gz";
#         subprocess.call(index, shell=True);


def clean_files(filename):
    """Check if intermediate but not final files exist and clean."""
    vcf_exists = os.path.exists(filename)
    gz_vcf_exists = os.path.exists(filename + '.gz')
    tbx_exists = os.path.exists(filename + '.gz.tbi')
    if vcf_exists and not tbx_exists:
        print('Deleting and restarting this file: ' + filename)
        os.remove(filename)
    if gz_vcf_exists and not tbx_exists:
        print('Deleting and restarting this file: ' + filename + '.gz')
        os.remove(filename + '.gz')


def split_compress_index(chrom, kiddo, fam_id):
    """Split trio VCF by chromosome and compress."""
    input_vcf = kiddo + '_trio.vcf.gz'
    filename = '{}/Illumina_WGS_{}_{}.vcf'.format(fam_id, fam_id, chrom)
    # check if already done
    if os.path.exists(filename + '.gz.tbi'):
        # check if empty file
        if os.stat(filename + '.gz.tbi').st_size == 0:
            print('deleting indexed ' + filename)
            # os.remove(filename + '.gz.tbi')
        else:
            print('Not re-running {} {}'.format(fam_id, chrom))
            return 'not_rerun_' + str(chrom)
    print('Starting {} {}'.format(fam_id, chrom))
    # delete intermediate/unfinished files
    clean_files(filename)
    # run actual commands
    split = 'time tabix -h {} {} > {}'.format(input_vcf, chrom, filename)
    sp.call(split, shell=True)
    compress = 'time bgzip ' + filename
    sp.call(compress, shell=True)
    index = 'time tabix -p vcf ' + filename + '.gz'
    sp.call(index, shell=True)
    return chrom


if __name__ == '__main__':
    trio_df = get_trio_df()
    # loop over kids here
    n_fams = trio_df.shape[0]
    for kiddo_ct in range(0, n_fams):
        kiddo = trio_df.loc[kiddo_ct, 'Child']
        fam_id = trio_df.loc[kiddo_ct, 'Fam_ID']
        print("Current ID {} from {}".format(kiddo, fam_id))
        # skip this trio if it hasn't been extracted from the main VCF
        if not os.path.exists(kiddo + '_trio.vcf.gz.tbi'):
            continue
        # create a directory per family ID if it doesn't exist
        if not os.path.exists(fam_id):
            print("Creating {} folder".format(fam_id))
            os.makedirs(fam_id)
        # get all chromosomes aka contigs from VCF
        tbx_handle = pysam.TabixFile(kiddo + '_trio.vcf.gz')
        contigs = tbx_handle.contigs
        # delete intermediate files that take up a lot of memory
        del tbx_handle
        gc.collect()
        split_compress_index_partial = partial(
            split_compress_index, kiddo=kiddo, fam_id=fam_id)
        chr_done = []
        for chrom in contigs:
            chr_done.append(split_compress_index_partial(chrom))
        # pool = mp.Pool(processes=5)
        print(chr_done)
        # chr_done = pool.map(split_compress_index_partial, contigs)
        # not_rerun_chr_list = [i for i in chr_done if 'not_rerun_' in i]
        # if len(not_rerun_chr_list) == 24:
        #     mv_cmd = 'mv {}_trio.vcf.g* done_2018_09_26/'.format(kiddo)
        #     print(mv_cmd)
        #     sp.call(mv_cmd, shell=True)
        # free up memory pool,
        del kiddo, fam_id, contigs, chr_done
        gc.collect()
