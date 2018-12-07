#!/usr/bin/env python3
# -*- coding: utf-8 -*-

r"""Classify variants as being in CpG sites or not using trinucleotide context.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-12-06
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

import pybedtools


def mod_start_end(feature):
    """Change start and end to trinucleotide space."""
    feature['start'] -= 1
    feature['stop'] += 1
    return feature


def write_trinuc(dnv_bed, out_loc):
    """Write the trinucleotide sequence context to a file."""
    with open(dnv_bed.seqfn) as seq_f, open(out_loc, 'w') as out_f:
        for line in seq_f:
            loc = line[1:].strip()
            if not loc.startswith('chr'):
                print(line)
                raise TypeError
            trinuc_line = next(seq_f).strip()
            outline = loc + '\t' + trinuc_line + '\n'
            _ = out_f.write(outline)
    return out_loc
    print(_)


def get_trinuc(bed_i, fa_loc):
    """Get trinucleotide context per DNV."""
    trinuc_out_loc = 'trinuc_out/{}_trinuc.txt'.format(bed_i[:-4])
    if os.path.exists(trinuc_out_loc):
        return 'trinuc already made for ' + bed_i
    dnv_bed = pybedtools.BedTool(bed_i)
    # change start column to start-1, end to end+1
    dnv_bed = dnv_bed.each(mod_start_end)
    dnv_bed = dnv_bed.saveas()
    # run getfasta using pybedtools wrapper
    fasta = pybedtools.example_filename(fa_loc)
    dnv_bed = dnv_bed.sequence(fi=fasta)
    # save output sequence
    write_trinuc(dnv_bed, trinuc_out_loc)
    return bed_i + ' trinuc done'


def get_cpg_isle(bed_i, cpg_island_bed):
    """Overlap DNVs with CpG islands."""
    cpg_isl_out_loc = 'CpG_islands/{}_cpg_isle.bed'.format(bed_i[:-4])
    if os.path.exists(cpg_isl_out_loc):
        return 'CpG islands already found for ' + bed_i
    dnv_bed = pybedtools.BedTool(bed_i)
    dnv_bed.intersect(cpg_island_bed).saveas(cpg_isl_out_loc)
    return bed_i + ' CpG island done'


pybedtools.set_tempdir('/sc/orga/projects/chdiTrios/Felix/tmp_pybedtools/')
os.chdir('/hpc/users/richtf01/longreadclustersequencing/data/gmkf2')
fa_loc = '/sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa'
cpg_island_bed = ('/hpc/users/richtf01/longreadclustersequencing/' +
                  'phasing_analysis/CpG_islands.bed')
bed_iter = glob.iglob('1-*_dnv.bed')
for bed_i in bed_iter:
    print(bed_i)
    get_trinuc(bed_i, fa_loc)
    get_cpg_isle(bed_i, cpg_island_bed)


#
