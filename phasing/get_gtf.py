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
IlluminaWhatshapVCFs/
python3 ~/longreadclustersequencing/phasing/get_gtf.py
Figure out how to get all files in sub-directories with glob

"""

import subprocess as sp
import multiprocessing as mp
import glob

# Script to get gtf files for all chromosomes for all IDs using multiprocessing

patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389",
             "1-04460", "1-04537", "1-05673", "1-05846"]


def get_gtf(ID):
    """Run Whatshap's gtf function to get contiguously phased variants."""
    for i in range(1, 23):
        command = ("whatshap stats --gtf=" + ID + "_chr" + str(i) +
                   "_phased.gtf " + ID + "_no_indels/" + ID + "_chr"
                   + str(i) + "_phased.vcf")
        sp.call(command, shell=True)
        print("Created gtf for " + ID + " chromosome " + str(i))


def get_ilmn_vcf_list():
    """Run Whatshap gtf for contiguously phased variants in Illumina data."""
    id_file_list = []
    for batch_i in ['1', '2', '3']:
        batch_f = 'Batch{}/*'.format(batch_i)
        id_file_list.extend([i for i in glob.iglob(batch_f)])
    ilmn_vcf_list = []
    for id_folder in id_file_list:
        ilmn_vcf_list.extend([i for i in glob.iglob(id_folder + '/*vcf')])
    return ilmn_vcf_list


def get_gtf_ilmn(vcf):
    """Run Whatshap gtf for contiguously phased variants in Illumina data."""
    gtf_cmd = "whatshap stats --gtf={}.gtf {}".format(vcf[:-4], vcf)
    print(gtf_cmd)
    # sp.call(gtf_cmd, shell=True)
    print("Created gtf for " + vcf)


if __name__ == '__main__':
    pool = mp.Pool(processes=5)
    # pool.map(get_gtf, patientID)
    ilmn_vcf_list = get_ilmn_vcf_list()
    print(ilmn_vcf_list)
    print(len(ilmn_vcf_list))
    pool.map(get_gtf_ilmn, ilmn_vcf_list)
