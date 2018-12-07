"""Run parent assignment algorithm on all patient IDs using phased VCFs.

Whatshap software (without indels flag)

module purge
module load python/3.5.0 py_packages/3.5
cd /hpc/users/richtf01/longreadclustersequencing/phasing
# python3 get_ID_dataframes.py --batch 1
python3

"""

# import multiprocessing as mp
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
    print('Number of patients overall: {}'.format(len(patientIDs)))
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
    print('Number of patients in batch: {}'.format(len(patientIDs)))
    # patientIDs.remove('1-05679')
    # if '1-04891' in patientIDs:
    #     patientIDs.remove('1-04891')
    #     # KeyError: 175492 on line 244: hap = curr_vcf[self.id][l_discon]
    return patientIDs


def write_missing_data(missing_list, missing_f):
    """Write list of missing files to `missing_f`."""
    if len(missing_list) > 0:
        with open(missing_f, 'a') as f:
            for i in missing_list:
                _ = f.write(i + '\n')
        print(_)


def write_problem_data(batch_ct, problem_pts):
    """Write list of problematic patients to a file."""
    key_error_f = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                   'PacbioProject/IlluminaWhatshapVCFs/keyerror_ids_b{}.txt')
    # Need to decide if writing or appending
    with open(key_error_f.format(batch_ct), 'w') as f:
        for id in problem_pts:
            _ = f.write(id + '\n')
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
    if os.path.exists(patient.out_f.format('')):
        print('already_done_' + patient.id)
        return 'already_done_' + patient.id
    if os.path.exists(patient.out_f.format('_incomplete')):
        print('previously_incomplete_' + patient.id)
        return 'previously_incomplete_' + patient.id
    home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                'PacbioProject/IlluminaWhatshapVCFs/')
    whatshap_prefix = (home_dir + 'Batch' + batch_i +
                       '/{}/{}_chr{}_phased')
    patient.illumina(whatshap_prefix)
    write_missing_data(
        patient.vcfs_todo, '{}vcfs_todo_b{}.txt'.format(home_dir, batch_i))
    write_missing_data(
        patient.gtfs_todo, '{}gtfs_todo_b{}.txt'.format(home_dir, batch_i))
    return patient.id


if __name__ == '__main__':
    # pool = mp.Pool(processes=5)
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
        get_illumina_GMKF2_dataframes, batch_i=batch_ct)
    done_pts = []
    problem_pts = []
    for ptID in patientIDs:
        print(ptID)
        try:
            done_pts.append(get_illumina_GMKF2_dataframes_partial(ptID))
        except KeyError:
            print('KeyError occurred with ' + ptID)
            problem_pts.append(ptID)
        except ValueError:
            print('ValueError occurred with ' + ptID)
            problem_pts.append(ptID)
    print('{} patients with KeyError or ValueError'.format(len(problem_pts)))
    print('{} patients completed successfully'.format(len(done_pts)))
    write_problem_data(batch_ct, problem_pts)
    # done_pts = pool.map(get_illumina_GMKF2_dataframes_partial, patientIDs)


"""Testing: figure out when and why the KeyErrors and ValueErrors occur

batch_ct = '2'
ptID = '1-06269'
get_illumina_GMKF2_dataframes(ptID, batch_ct)

patientIDs = get_patient_ids(batch_ct)
get_illumina_GMKF2_dataframes_partial = partial(
    get_illumina_GMKF2_dataframes, batch_i=batch_ct)
done_pts = []
problem_pts = []
for ptID in patientIDs:
    print(ptID)
    try:
        done_pts.append(get_illumina_GMKF2_dataframes_partial(ptID))
    except KeyError:
        print('KeyError occured with ' + ptID)
        problem_pts.append(ptID)


key_error_f = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
               'PacbioProject/IlluminaWhatshapVCFs/keyerror_ids_b{}.txt')
with open(key_error_f.format(batch_ct), 'w') as f:
    for id in problem_pts:
        _ = f.write(id + '\n')


batch_i = str(1)
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

# confirm all IDs are accounted for. 1 missing, that's fine
# 123 + 194 + 38, 356


/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/IlluminaWhatshapVCFs/Batch2/CG0017-8871/
1-12147_chr1_phased.vcf

chr1:8867755

"""
