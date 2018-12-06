#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run sorting-hat on all outputs.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-12-05
:Copyright: 2018, Allison Seiden, Felix Richter
:License: CC BY-SA

cd /hpc/users/richtf01/longreadclustersequencing/data/gmkf2

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix/0.2.6 bedtools/2.27.1
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
python3

"""

import os
import glob
import re
import subprocess as sp


def run_sort_hat(bed_i, fa_loc, rmsk):
    """Run sorting-hat for a single bed file."""
    cmd = 'time sorting_hat --bed {} --fasta {} --repeat {} --output {}'
    out_f = 'sorting_hat_out/' + re.sub('bed$', 'txt', bed_i)
    if os.path.exists(out_f):
        print(out_f + ' already exists')
        return out_f + ' already exists'
    cmd = cmd.format(bed_i, fa_loc, rmsk, out_f)
    print(cmd)
    sp.call(cmd, shell=True)
    return bed_i


os.chdir('/hpc/users/richtf01/longreadclustersequencing/data/gmkf2')

bed_iter = glob.iglob('1-*_dnv.bed')
fa_loc = '/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa'
rmsk = '/hpc/users/richtf01/longreadclustersequencing/data/repeat_masker.txt'

bed_done_list = []
for bed_i in bed_iter:
    print(bed_i)
    bed_done = run_sort_hat(bed_i, fa_loc, rmsk)
    bed_done_list.append(bed_done)


#
#
