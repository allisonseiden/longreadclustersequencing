"""Run parent assignment algorithm on all patient IDs using phased VCFs.

Whatshap software (without indels flag)

module purge
module load python/3.5.0 py_packages/3.5
cd /hpc/users/richtf01/longreadclustersequencing/phasing
# python3 get_ID_dataframes.py
python3

"""

import multiprocessing as mp
import re
import os
from functools import partial
import argparse


from PhasedData import PhasedData
from utils import get_trio_df, get_done_files, get_batch_pt_ids


def get_patient_ids(batch_i):
    """Get list of patients IDs for DNV phasing."""
    # get previously completed IDs
    done_list = get_done_files()
    patientIDs = list(set([re.sub('.*_fullphased/|_chr.*', '', i)
                           for i in done_list[1:]]))
    print(len(patientIDs))
    # patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389',
    #               '1-04460', '1-04537', '1-05443', '1-05673', '1-05846']
    # patientIDs = ['1-06149', '1-05794', '1-05935', '1-05860', '1-05423']
    # only keep IDs in a specific batch
    b_id_list = get_batch_pt_ids(batch_i)
    b_fam_id_list = []
    trio_df = get_trio_df()
    for b_id in b_id_list:
        b_fam_id_list.append(
            trio_df.loc[trio_df.Child == b_id][
                'Fam_ID'].to_string(index=False))
    patientIDs = [i for i in patientIDs if i in b_fam_id_list]
    print(len(patientIDs))
    # patientIDs.remove('1-05679')
    patientIDs.remove('1-04891')
    # KeyError: 175492 on line 244: hap = curr_vcf[self.id][l_discon]
    return patientIDs


def write_missing_data(missing_list, missing_f):
    """Write list of missing files to `missing_f`."""
    if len(missing_list) > 0:
        with open(missing_f, 'a') as f:
            for i in missing_list:
                _ = f.write(i + '\n')
        print(_)


def get_pacbio_dataframes(ID):
    """Get phased de novo variants for PacBio data."""
    patient = PhasedData(ID)
    patient.pacbio()


def get_illumina_dataframes(ID):
    """Get phased de novo variants for Illumina data."""
    patient = PhasedData(ID)
    whatshap_prefix = ('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' +
                       '{}_illumina/{}_chr{}_phased')
    patient.illumina(whatshap_prefix)


def get_illumina_GMKF2_dataframes(ID, batch_i):
    """Get phased de novo variants for Illumina data."""
    # provide trio_df if VCF IDs are not the family IDs
    trio_df = get_trio_df()
    patient = PhasedData(ID, trio_df, home_dir='/hpc/users/richtf01/')
    if os.path.exists('phased_data/' + patient.id + '_dataframe.txt'):
        return 'already_done_' + patient.id
    home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                'PacbioProject/IlluminaWhatshapVCFs/')
    whatshap_prefix = (home_dir + 'Batch' + batch_i +
                       '/{}/{}_chr{}_phased')
    patient.illumina(whatshap_prefix)
    write_missing_data(
        patient.vcfs_todo, '{}vcfs_todo_b{}.txt'.format(home_dir, batch_i))
    write_missing_data(
        patient.gtfs_todo, '{}gtfs_todo{}.txt'.format(home_dir, batch_i))
    return patient.id


if __name__ == '__main__':
    pool = mp.Pool(processes=5)
    parser = argparse.ArgumentParser()
    parser.add_argument('--batch', default='1', type=str,
                        choices=['1', '2', '3'], help='Pick the batch')
    args = parser.parse_args()
    batch_ct = args.batch
    # comment out the one which you're not using at the moment, only call one
    # of these at a time
    # pool.map(get_pacbio_dataframes, patientIDs)
    # pool.map(get_illumina_dataframes, patientIDs)
    patientIDs = get_patient_ids(batch_ct)
    get_illumina_GMKF2_dataframes_partial = partial(
        get_illumina_GMKF2_dataframes, batch_ct)
    pool.map(get_illumina_GMKF2_dataframes, patientIDs)


"""Testing

batch_i = str(2)
patientIDs = get_patient_ids(batch_i)
patientIDs[0]
get_illumina_GMKF2_dataframes(patientIDs[0], batch_i)

trio_df = get_trio_df()
ID = '1-05794'
patient = PhasedData(ID, trio_df, home_dir='/hpc/users/richtf01/')
whatshap_prefix = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                   'PacbioProject/IlluminaWhatshapVCFs/Batch' + batch_i +
                   '/{}/{}_chr{}_phased')

patient.illumina(whatshap_prefix)

pool = mp.Pool(processes=3)
pool.map(get_illumina_GMKF2_dataframes, patientIDs)

"""
