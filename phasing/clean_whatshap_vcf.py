#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Clean whatshap output VCF.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-10-15
:Copyright: 2018, Allison Seiden
:License: CC BY-SA


module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/
cd ~/longreadclustersequencing/phasing
python3

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/
python3 ~/longreadclustersequencing/phasing/clean_whatshap_vcf.py --batch 1

"""

import glob
import re
import os
import subprocess as sp
import argparse
# import multiprocessing as mp

from utils import get_done_files

# read in phased_vcf_filtered_done.txt
# skip those IDs that are done
# move the original whatshap VCF to a 'fullVCF' folder
# remove 0/0 and ./. lines from whatshap VCF
# remove if ALL 3 are 0/0 or ./., or only the proband?

vcf_cols = ['CHROM', 'POS', 'ID', 'REF',
            'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
            'ID', 'MOM', 'DAD']

chrom_len_dict = {'chr1': 5519702, 'chr2': 5320494, 'chr3': 4315874,
                  'chr4': 4201162, 'chr5': 3890615, 'chr6': 3677158,
                  'chr7': 3624114, 'chr8': 3295633, 'chr9': 2925869,
                  'chr10': 3124593, 'chr11': 3140594, 'chr12': 3032469,
                  'chr13': 2468214, 'chr14': 1951375, 'chr15': 1881407,
                  'chr16': 2064351, 'chr17': 1902520, 'chr18': 1851434,
                  'chr19': 1449836, 'chr20': 1568104,
                  'chr21': 991855, 'chr22': 1038052, 'chrX': 1,
                  'chrY': 1}


done_list = get_done_files()


def check_null(format_col):
    """Check for removal evidence."""
    return format_col.startswith('0/0') or format_col.startswith('./.')


def clean_vcf(vcf_loc, clean_vcf):
    """Write clean VCF to new, hopefully much smaller file."""
    with open(vcf_loc, 'r') as f, open(clean_vcf, 'w') as out_f:
        line = next(f)
        while line.startswith('##'):
            _ = out_f.write(line)
            line = next(f)
        # write the header line
        _ = out_f.write(line)
        count = 0
        for line in f:
            line_dict = dict(zip(vcf_cols, line.strip().split('\t')))
            rm_id = check_null(line_dict['ID'])
            """DECIDE IF PROBAND ONLY OR ALL."""
            rm_mom = check_null(line_dict['MOM'])
            rm_dad = check_null(line_dict['DAD'])
            # if rm_id:
            count += 1
            if (count % 500000) == 0:
                print("Fraction done: {}".format(
                    round(count/chrom_len_dict[line_dict['CHROM']], 3)))
            if rm_id and rm_mom and rm_dad:
                continue
            else:
                _ = out_f.write(line)
        print(count)
        frac_done = round(count/chrom_len_dict[line_dict['CHROM']], 3)
        print("Final raction done: {}".format(frac_done))
        if frac_done != 1.000:
            print(vcf_loc)
            print(frac_done)
            raise ValueError('Non-1 fraction of chr done for' + vcf_loc)
    return vcf_loc
    print(_)  # only have this to get rid of pylint warning


def get_ilmn_vcf_list(batch_list):
    """Get a list of illumina VCFs that have been phased."""
    id_file_list = []
    # '1',, '2', '3'
    for batch_i in batch_list:
        batch_f = 'Batch{}/*'.format(batch_i)
        id_file_list.extend([i for i in glob.iglob(batch_f)])
    ilmn_vcf_list = []
    for id_folder in id_file_list:
        ilmn_vcf_list.extend([i for i in glob.iglob(id_folder + '/*vcf')])
    return ilmn_vcf_list


def process_vcf(vcf):
    """Move and clean VCF."""
    final_vcf_loc = vcf
    done_list = get_done_files()
    if not os.path.exists(final_vcf_loc):
        return None
    full_vcf_loc = re.sub('Batch[1-3]/.*/',
                          '../IlluminaWhatshapVCFs_fullphased/', vcf)
    if not (full_vcf_loc in done_list):
        # move file to phased VCF directory if it does not already exist
        if not os.path.exists(full_vcf_loc):
            mv_cmd = 'mv {} {}'.format(final_vcf_loc, full_vcf_loc)
            print(mv_cmd)
            sp.call(mv_cmd, shell=True)
        vcf_done = clean_vcf(full_vcf_loc, clean_vcf=final_vcf_loc)
        print(vcf_done)
        # After running append to phased_vcf_filtered_done.txt
        # Should only append if process ended successfully
        log_f_loc = 'vcf_cleaning_tracker.txt'
        with open(log_f_loc, 'a') as log_f:
            log_f.write(vcf_done + '\n')
        print('removing', full_vcf_loc)
        os.remove(full_vcf_loc)
        return full_vcf_loc
    else:
        print('already completed', full_vcf_loc)
        return 'done: ' + full_vcf_loc


# def clean_full_phased_dir(vcf):
#     """Delete full VCF."""
#     final_vcf_loc = vcf
#     done_list = get_done_files()
#     if not os.path.exists(final_vcf_loc):
#         return None
#     full_vcf_loc = re.sub('Batch[1-3]/.*/',
#                           '../IlluminaWhatshapVCFs_fullphased/', vcf)
#     if (full_vcf_loc in done_list) and os.path.exists(full_vcf_loc):
#         # if os.path.exists(full_vcf_loc):
#         os.remove(full_vcf_loc)
#         print('removing', full_vcf_loc)


if __name__ == '__main__':
    os.chdir('/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
             'IlluminaWhatshapVCFs/')
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch', default='1', type=str,
                        choices=['1', '2', '3'], help='Pick the batch')
    args = parser.parse_args()
    batch_ct = args.batch
    ilmn_vcf_list = get_ilmn_vcf_list(batch_ct)
    done_list = get_done_files()
    print(len(ilmn_vcf_list))
    print(len(done_list))
    # pool = mp.Pool(processes=5)
    # clean_list = pool.map(process_vcf, ilmn_vcf_list)
    for ilmn_vcf in ilmn_vcf_list:
        print(ilmn_vcf)
        process_vcf(ilmn_vcf)


"""
os.chdir('/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
         'IlluminaWhatshapVCFs/')
batch_list = ['1', '2', '3'][0]
ilmn_vcf_list = get_ilmn_vcf_list(batch_list)
done_list = get_done_files()
print(len(ilmn_vcf_list))
len(done_list)

for ilmn_vcf in ilmn_vcf_list:
    print(ilmn_vcf)
    process_vcf(ilmn_vcf)

pool = mp.Pool(processes=5)
clean_list = pool.map(process_vcf, ilmn_vcf_list)

"""

"""

# figure out what's going on with:
rm IlluminaWhatshapVCFs/Batch2/CG0016-3456/1-07336_chr4_phased.vcf
rm IlluminaWhatshapVCFs_fullphased/1-07336_chr4_phased.vcf
vim IlluminaWhatshapVCFs/vcf_cleaning_tracker.txt
DOES NOT ACTUALLY LOOK LIKE THERE'S A MISTAKE...

# fixing a mistake
from utils import get_batch_pt_ids, get_trio_df, get_done_files

batch_ct = '2'
patientID_list = get_batch_pt_ids(batch_ct)
trio_df = get_trio_df()
re_str = '../IlluminaWhatshapVCFs_fullphased/|_phased.vcf'

# '../IlluminaWhatshapVCFs_fullphased/1-02284_chr1_phased.vcf'
# 'Batch3/CG0005-3923/1-02468_chr14_phased.vcf'

count = 0
for full_vcf_loc in done_list:
    # get final_vcf_loc from full_vcf_loc somehow
    f_info_list = re.sub(re_str, '', full_vcf_loc).split('_')
    sample_id = [trio_df[trio_df.Fam_ID == f_info_list[0]]
                 ['Child'].to_string(index=False)][0]
    if not (sample_id in patientID_list):
        continue
    final_vcf_loc = 'Batch{}/{}/{}_{}_phased.vcf'.format(
        batch_ct, sample_id, f_info_list[0], f_info_list[1])
    if os.path.exists(final_vcf_loc):
        # print(final_vcf_loc)
        count += 1
    if not os.path.exists(final_vcf_loc):
        mv_cmd = 'mv {} {}'.format(full_vcf_loc, final_vcf_loc)
        print(mv_cmd)
        # sp.call(mv_cmd, shell=True)

for ilmn_vcf in ilmn_vcf_list:
    print(ilmn_vcf)
    # clean_full_phased_dir(ilmn_vcf)

"""

"""OTHER SANITY CHECKS
# All batch3 and batch2 files in done_files should have their non-zero VCF
# in the Batch2/ and Batch3/ directories

from utils import get_batch_pt_ids, get_trio_df, get_done_files

os.chdir('/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/' +
         'IlluminaWhatshapVCFs/')
batch_ct = '2'
done_list = get_done_files()
patientID_list = get_batch_pt_ids(batch_ct)
trio_df = get_trio_df()
re_str = '../IlluminaWhatshapVCFs_fullphased/|_phased.vcf'
correct_loc_ct = 0
done_in_batch_ct = 0
log_f_loc = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
             'PacbioProject/IlluminaWhatshapVCFs/vcf_cleaning_tracker_b{}.txt')

with open(log_f_loc.format(batch_ct), 'w') as f:
    for full_vcf_loc in done_list:
        # get final_vcf_loc from full_vcf_loc
        f_info_list = re.sub(re_str, '', full_vcf_loc).split('_')
        sample_id = [trio_df[trio_df.Fam_ID == f_info_list[0]]
                     ['Child'].to_string(index=False)][0]
        if not (sample_id in patientID_list):
            continue
        final_vcf_loc = 'Batch{}/{}/{}_{}_phased.vcf'.format(
            batch_ct, sample_id, f_info_list[0], f_info_list[1])
        done_in_batch_ct += 1
        # print(final_vcf_loc)
        if os.path.exists(final_vcf_loc):
            correct_loc_ct += 1
            _ = f.write(full_vcf_loc + '\n')
        else:
            tmp_loc = 'tmp_to_keep/Batch{}_{}_{}_{}_phased.vcf'.format(
                batch_ct, sample_id, f_info_list[0], f_info_list[1])
            mv_cmd = 'mv {} {}'.format(full_vcf_loc, tmp_loc)
            print(mv_cmd)
            sp.call(mv_cmd, shell=True)


print(done_in_batch_ct, correct_loc_ct)

cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/\
IlluminaWhatshapVCFs/
mv ../IlluminaWhatshapVCFs_fullphased/1-04961_chr7_phased.vcf \
Batch2/CG0013-8579/1-04961_chr7_phased.vcf
mv ../IlluminaWhatshapVCFs_fullphased/1-04961_chr8_phased.vcf \
Batch2/CG0013-8579/1-04961_chr8_phased.vcf
# delete these two files from done_files

mv vcf_cleaning_tracker.txt vcf_cleaning_tracker_retired_11_04.txt
cat vcf_cleaning_tracker_b* > vcf_cleaning_tracker.txt

# All batch3 and batch2 files NOT in done_files should be the original size

# What is the deal with the extra Tb of IlluminaWhatshapVCFs_fullphased files?
"""
