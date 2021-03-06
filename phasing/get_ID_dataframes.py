"""Run parent assignment algorithm on all patient IDs using phased VCFs.

Whatshap software (without indels flag)

module purge
module load python/3.5.0 py_packages/3.5
cd /hpc/users/richtf01/longreadclustersequencing/phasing
# python3 get_ID_dataframes.py --batch 3
python3

"""

# import multiprocessing as mp
import os
from functools import partial
import argparse


from PhasedData import PhasedData
from utils import get_trio_df, get_patient_ids


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
                       '{}_illumina_2019/{}_chr{}_phased')
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
        # try:
        done_pts.append(get_illumina_GMKF2_dataframes_partial(ptID))
        # except KeyError:
        #     print('KeyError occurred with ' + ptID)
        #     problem_pts.append(ptID)
        # except ValueError:
        #     print('ValueError occurred with ' + ptID)
        #     problem_pts.append(ptID)
    # print('{} patients with KeyError or ValueError'.format(len(problem_pts)))
    # write_problem_data(batch_ct, problem_pts)
    print('{} patients completed successfully'.format(len(done_pts)))
    # done_pts = pool.map(get_illumina_GMKF2_dataframes_partial, patientIDs)


"""Testing for Ilmn 10:
ID = '1-01019'
patient = PhasedData(ID)
whatshap_prefix = ('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' +
                   '{}_illumina_2019/{}_chr{}_phased')
patient.illumina(whatshap_prefix)

## looping (ignoring '1-01019', '1-00801', )
ID_list = ['1-03897', '1-04190', '1-04389',
           '1-04460', '1-04537', '1-05443', '1-05673',
           '1-05846']

for ID_i in ID_list:
    get_illumina_dataframes(ID_i)

"""

"""Testing: figure out when and why the KeyErrors and ValueErrors occur

batch_ct = '2'
ptID = '1-06269'
get_illumina_GMKF2_dataframes(ptID, batch_ct)

## Troubleshooting specific steps
trio_df = get_trio_df()
patient = PhasedData(ptID, trio_df, home_dir='/hpc/users/richtf01/')
home_dir = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
            'PacbioProject/IlluminaWhatshapVCFs/')
whatshap_prefix = (home_dir + 'Batch' + batch_ct +
                   '/{}/{}_chr{}_phased')
patient.create_vcf_dictionary(whatshap_prefix)
# clean_gtf_output is only for Illumina GMKF (to check if done)
patient.clean_gtf_output(whatshap_prefix)
patient.create_dnvs_dictionary()


# figure out what's going on at chr3 where dnv == 174059053
chr_bounds = {}
# Collect correct phased VCF file for the current chromosome
chromosome = 'chr3'
curr_vcf = patient.vcf_dfs[chromosome]
# Collect start and end positions of haplotype blocks from GTF file for
# current chromosome
start_list = patient.gtf_dfs[chromosome]['Start'].tolist()
end_list = patient.gtf_dfs[chromosome]['End'].tolist()
end_list[-1]
dnv = 174059053

patient.fill_bounds_dictionary()
patient.find_variants_for_phasing(7)
patient.assign_to_parent()
patient.convert_to_dataframe()

# confirm all IDs are accounted for. 1 missing, that's fine
# 123 + 194 + 38, 356

"""
