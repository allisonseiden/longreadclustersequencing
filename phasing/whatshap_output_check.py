#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Check whatshap output length, then filter for only lines with vars in ID

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-19
:Copyright: 2018, Allison Seiden
:License: CC BY-SA


module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch2/
python3 ~/longreadclustersequencing/phasing/whatshap_output_check.py --batch 2

"""

import argparse
import subprocess as sp

from utils import get_batch_pt_ids, get_trio_df

# def mapcount(filename):
#     """Source: https://stackoverflow.com/a/850962"""
#     f = open(filename, "r+")
#     buf = mmap.mmap(f.fileno(), 0)
#     lines = 0
#     readline = buf.readline
#     while readline():
#         lines += 1
#     return lines


def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)


def rawgencount(filename):
    """Source: https://stackoverflow.com/a/27518377 ."""
    f = open(filename, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum(buf.count(b'\n') for buf in f_gen)


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
    patientID_list = get_batch_pt_ids(batch_ct)
    trio_df = get_trio_df()
    chr_dict = {}
    count = 0
    for i in range(1, 23):
        len_dict = {}
        chr_dict[i] = len_dict
        for ID in patientID_list:
            fam_id = trio_df.loc[
                trio_df.Child == ID]['Fam_ID'].to_string(index=False)
            phase_f = '{}/{}_chr{}_phased.vcf'.format(ID, fam_id, i)
            len_dict[phase_f] = rawgencount(phase_f)
            count += 1
            if count > 10:
                break
    print(chr_dict)

#