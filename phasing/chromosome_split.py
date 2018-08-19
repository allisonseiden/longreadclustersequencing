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

Once done, check that last line in every file is the same

"""

import subprocess
import os
import multiprocessing as mp
from functools import partial

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


def split_compress_index(chrom, kiddo, fam_id):
    """Split trio VCF by chromosome and compress."""
    print(kiddo)
    print(fam_id)
    input_vcf = kiddo + '_trio.vcf.gz'
    filename = '{}/Illumina_WGS_{}_chr{}.vcf'.format(fam_id, fam_id, chrom)
    print(input_vcf)
    print(filename)
    split = 'tabix -h {} chr{} > {}'.format(input_vcf, chrom, filename)
    print(split)
    # subprocess.call(split, shell=True)
    # compress = 'bgzip ' + filename
    # subprocess.call(compress, shell=True)
    # index = 'tabix -p vcf ' + filename + '.gz'
    # subprocess.call(index, shell=True)


if __name__ == '__main__':
    trio_df = get_trio_df()
    # loop over kids here
    kiddo = trio_df.loc[0, 'Child']
    fam_id = trio_df.loc[0, 'Fam_ID']
    # create a directory per family ID if it doesn't exist
    if not os.path.exists(fam_id):
        print("Creating", fam_id)
        os.makedirs(fam_id)
    split_compress_index_partial = partial(
        split_compress_index, kiddo=kiddo, fam_id=fam_id)
    pool = mp.Pool(processes=5)
    pool.map(split_compress_index_partial, range(1, 22))
