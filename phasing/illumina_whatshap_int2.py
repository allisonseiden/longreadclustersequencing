"""Run whatshap with initial 10 triosself.

:Authors: Allison Seiden, Felix Richter
:Date: 2019-01-31
:Copyright: 2018, Allison Seiden
:License: CC BY-SA

Description: Run Whatshap with indels flag using Illumina data
    10 patient IDs, uses multiprocessing.

cd ~
module purge
module load samtools/1.8 bcftools/1.7 tabix
module load python/3.5.0 py_packages/3.5
source venv_phasing/bin/activate
python3

"""

import subprocess as sp
import multiprocessing as mp
import os

# df = pd.read_table(
#     '/hpc/users/seidea02/longreadclustersequencing/data/' +
#     'Illumina_crams_batch2.txt', names=['ID'])
# patientID = df['ID'].tolist()

# 1-00801 was copied to an andy sharp directory I don't have access to
patientID = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389'][1:5]
# patientID = ['1-04460', '1-04537', '1-05443', '1-05673', '1-05846'][1:5]


def illumina_whatshap(ID):
    """Run illumina whatshap command for original 10 trios."""
    mkdir = 'mkdir -p ' + ID + '_illumina_2019'
    print(mkdir)
    print(os.getcwd())
    sp.call(mkdir, shell=True)
    for i in range(1, 23):
        vcf_filename = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                        'PacbioProject/DNV_calls/VCF/TrioVCF/' + ID + '/' +
                        ID + '_chr' + str(i) + '.vcf.gz')
        bam_filename = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                        'PacbioProject/IlluminaHg38Crams/' + ID +
                        '.hg38.dedup.clean.recal.cram')
        out_f = '{}_illumina_2019/{}_chr{}_phased.vcf'.format(ID, ID, i)
        command = 'whatshap phase --sample=' + ID + ' --ignore-read-groups'
        command += ' --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa'
        command += ' --indels -o ' + out_f
        command += ' ' + vcf_filename + ' ' + bam_filename
        print(command)
        sp.call(command, shell=True)
        print('======Sucessfully ran whatshap for ' + ID +
              ' on chromosome ' + str(i))


if __name__ == '__main__':
    home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                'PacbioProject/WhatshapVCFs/')
    os.chdir(home_dir)
    pool = mp.Pool(processes=5)
    pool.map(illumina_whatshap, patientID)


"""Testing:
home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
            'PacbioProject/WhatshapVCFs/')
os.chdir(home_dir)
os.getcwd()

illumina_whatshap(patientID[1])

home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
            'PacbioProject/WhatshapVCFs/')
os.chdir(home_dir)
os.getcwd()
pool = mp.Pool(processes=2)
pool.map(illumina_whatshap, patientID)

"""
