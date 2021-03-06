#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Check whatshap output length and delete if not consistent.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-19
:Copyright: 2018, Allison Seiden
:License: CC BY-SA

Pipeline overview:
 - whasthap_bsub.sh (or illumina_whatshap_int1.py)
 - whasthap_output_check.py
 - clean_whatshap_vcf.py
 - get_gtf.py
 - get_ID_dataframes.py


module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch2/
python3 ~/longreadclustersequencing/phasing/whatshap_output_check.py --batch 2

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch3/
python3 ~/longreadclustersequencing/phasing/whatshap_output_check.py --batch 3

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/Batch1/
python3 ~/longreadclustersequencing/phasing/whatshap_output_check.py --batch 1

cd ~/longreadclustersequencing/phasing/
python

"""

import argparse
import subprocess as sp
import os

from utils import get_batch_pt_ids, get_trio_df, get_done_files

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


def remove_incomplete_files(len_dict):
    """Remove files that are not complete."""
    # for chr_i, len_dict in chr_dict.items():
    max_len = max(len_dict.values())
    print('max length:')
    print(max_len)
    done_files = 0
    # always remove zero-length files
    if max_len == 0:
        max_len = 1
    for f, f_len in len_dict.items():
        if f_len != max_len:
            print('removing {} w length {}'.format(f, f_len))
            # os.remove(f)
        else:
            done_files += 1
    print('files kept: {}'.format(done_files))


def check_and_rm_files(patientID_list, trio_df, done_list):
    """Check file length and remove if not complete."""
    chr_dict = {}
    for i in range(1, 23):
        count = 0
        len_dict = {}
        chr_dict[i] = len_dict
        for ID in patientID_list:
            fam_id = trio_df.loc[
                trio_df.Child == ID]['Fam_ID'].to_string(index=False)
            phase_f = '/{}_chr{}_phased.vcf'.format(fam_id, i)
            if any([phase_f in i for i in done_list]):
                print('Already decreased size of ' + phase_f)
                continue
            phase_f = ID + phase_f
            if os.path.isfile(phase_f):
                print(phase_f)
                len_dict[phase_f] = rawgencount(phase_f)
                print(len_dict[phase_f])
                count += 1
                # if count % 5 == 0:
                #     print(count)
            # if count > 5:
            #     break
        print('total files: {}'.format(count))
        # only remove files if there are enough representative chromosomes
        # where you think the maximum length is truly the max
        if len(len_dict) > 0:
            remove_incomplete_files(len_dict)
        else:
            print(len_dict)


def clean_old_vcfs(patientID_list, trio_df):
    """Move old VCF to directory for archiving if phasing is done."""
    vcf_dir = ('/sc/orga/projects/chdiTrios/' +
               'WGS_Combined_2017/PacbioProject/' +
               'GMKF_TrioVCFs')
    for ID in patientID_list:
        fam_id = trio_df.loc[
            trio_df.Child == ID]['Fam_ID'].to_string(index=False)
        # split_chr_done_2018_09_26
        new_dir = '{}/split_chr_done_2018_11_03/{}/'.format(vcf_dir, fam_id)
        dir_cmd = 'mkdir -p ' + new_dir
        print(dir_cmd)
        sp.call(dir_cmd, shell=True)
        for i in range(1, 23):
            # """
            phase_f = '{}/{}_chr{}_phased.vcf'.format(ID, fam_id, i)
            # Moving the corresponding split chromosome file
            # to a 'done' location
            vcf_old = ('{}/{}/Illumina_WGS_{}' +
                       '_chr{}.vcf').format(
                       vcf_dir, fam_id, fam_id, i)
            mv_cmd = 'mv {}.g* {}'.format(vcf_old, new_dir)
            phase_exists = os.path.isfile(phase_f)
            vcf_old_exists = os.path.isfile(vcf_old + '.gz.tbi')
            if phase_exists and vcf_old_exists:
                print(mv_cmd)
                sp.call(mv_cmd, shell=True)
            """
            # reversing a mistake...
            vcf_new = ('{}/split_chr_done_2018_09_26/{}/Illumina_WGS_{}' +
                       '_chr{}.vcf').format(
                       vcf_dir, fam_id, fam_id, i)
            old_dir = '{}/{}/'.format(vcf_dir, fam_id)
            mv_cmd = 'mv {}.g* {}'.format(vcf_new, old_dir)
            vcf_old_exists = os.path.isfile(vcf_new + '.gz.tbi')
            if vcf_old_exists:
                print(mv_cmd)
                sp.call(mv_cmd, shell=True)
            # """


if __name__ == '__main__':
    # pool = mp.Pool(processes=2)
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch', default=1, type=int,
                        choices=[1, 2, 3], help='Pick the batch')
    args = parser.parse_args()
    batch_ct = args.batch
    os.chdir('/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
             'IlluminaWhatshapVCFs/Batch' + str(batch_ct))
    patientID_list = get_batch_pt_ids(batch_ct)
    trio_df = get_trio_df()
    done_list = get_done_files()
    # Batch2/CG0012-6043/1-05794 to test done_list
    check_and_rm_files(patientID_list, trio_df, done_list)
    clean_old_vcfs(patientID_list, trio_df)

"""Testing

batch_ct = 3
os.chdir('/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
         'IlluminaWhatshapVCFs/Batch' + str(batch_ct))
patientID_list = get_batch_pt_ids(batch_ct)
trio_df = get_trio_df()
done_list = get_done_files()
# patientID_list = [i for i in patientID_list if i is '1-05794']
# patientID_list = [trio_df[trio_df.Fam_ID == '1-05794']
#                   ['Child'].to_string(index=False)]
# Batch2/CG0012-6043/1-05794 to test done_list

## For pacbio
patientID_list = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389'][1:5]
home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
            'PacbioProject/WhatshapVCFs/')
os.chdir(home_dir)
check_and_rm_files(patientID_list, trio_df, done_list)
clean_old_vcfs(patientID_list, trio_df) # not working (wrong filenames/dirs)
"""
#
