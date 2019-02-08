#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Replication timing overlap with DNVs.

:Authors: Felix Richter
:Date: 2019-01-07
:Copyright: 2019, Felix Richter
:License: CC BY-SA

cd /sc/orga/projects/chdiTrios/Felix/repli_seq_data/
cd ~/longreadclustersequencing
module load ucsc-utils/2015-04-07
module load python/3.5.0 py_packages/3.5
python3

"""

import glob
import re
import subprocess


"""Downloading the ENCODE data:
# DL links from:
https://www.encodeproject.org/matrix/?type=Experiment&y.limit=&assay_title=Repli-seq&assay_slims=Replication+timing&biosample_ontology.cell_slims=stem+cell
# ignore bams and fastq.gz files
grep -v 'bam$\|fastq.gz$' files\(1\).txt > dl_files.txt
# run download command
xargs -L 1 curl -O -L < dl_files.txt
"""

# Define file names
home_dir = '/hpc/users/richtf01/longreadclustersequencing/'
in_f_loc = home_dir + 'data/dnvs_2019_02_07.bed'
in_f_loc = home_dir + 'literature/jonsson_decode_2017/indels.bed'

"""create unique line IDs for every file."""
lineid_f_loc = in_f_loc[:-4] + '_lineid.bed'
lineid_f_loc = re.sub('/data/', '/data/repliseq_anno/', lineid_f_loc)
lineid_f_loc = re.sub(
    'decode_2017/', 'decode_2017/repliseq_anno/', lineid_f_loc)
with open(in_f_loc, 'r') as in_f, open(lineid_f_loc, 'w') as out_f:
    for line in in_f:
        line_list = line.strip().split("\t")
        var_id = '.'.join(line_list)
        line_list = line_list[0:3]
        line_list.append(var_id)
        # print(line_list)
        out_line = "\t".join(line_list) + "\n"
        _ = out_f.write(out_line)

"""sort file."""
sorted_f = lineid_f_loc[:-4] + '_sorted.bed'
sort_cmd = 'time sort -V -k1,1 -k2,2 {} | uniq > {}'.format(
    lineid_f_loc, sorted_f)
print(sort_cmd)
subprocess.call(sort_cmd, shell=True)

"""lift over to hg19."""
bed_dir = "/hpc/users/richtf01/chdiTrios/Felix/wgs/bed_annotations/"
lift_f = sorted_f[:-4] + '_hg19.bed'
lift_f_unmapped = lift_f[:-4] + '_unmapped.bed'
# hg38ToHg19 hg19ToHg38
chain_f = bed_dir + "hg38ToHg19.over.chain"
liftOverCommand = "time liftOver -multiple -minMatch=0.1 {} {} {} {}".format(
    sorted_f, chain_f, lift_f, lift_f_unmapped)
print(liftOverCommand)
subprocess.call(liftOverCommand, shell=True)

"""Overlap with replication timing data."""

bw_f_glob = '/sc/orga/projects/chdiTrios/Felix/repli_seq_data/*.bigWig'
bw_f_loc_list = [i for i in glob.iglob(bw_f_glob)]
bw_f_loc_i = bw_f_loc_list[0]

for bw_f_loc_i in bw_f_loc_list:
    bw_f_name_i = re.sub('.*/|.bigWig', '', bw_f_loc_i)
    bw_cmd = ('time bigWigAverageOverBed {} {} {}_{}.tab ' +
              '-bedOut={}_{}.bed').format(
        bw_f_loc_i, lift_f, lift_f[:-4], bw_f_name_i,
        lift_f[:-4], bw_f_name_i)
    print(bw_cmd)
    subprocess.call(bw_cmd, shell=True)


"""End of script."""
