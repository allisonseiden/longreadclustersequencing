#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions called in multiple files.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-19
:Copyright: 2018, Allison Seiden
:License: CC BY-SA

"""

import re
import pandas as pd


def get_trio_df():
    """Load Illumina trio dataframe."""
    trio_df = pd.read_table('/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_' +
                            '2017/GMKF_Seidman_CongenitalHeartDisease_WGS.ped',
                            names=['Fam_ID', 'Child', 'Father', 'Mother'],
                            sep='\t', comment='#')
    return trio_df


def get_batch_pt_ids(batch_ct):
    """Load data on the correct batch."""
    df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/' +
                       'Illumina_crams_batch{}.txt'.format(batch_ct),
                       names=['ID'])
    patientID = df['ID'].tolist()
    return patientID


def get_done_files():
    """Get list of files completed."""
    done_list = []
    log_f_loc = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                 'PacbioProject/IlluminaWhatshapVCFs/vcf_cleaning_tracker.txt')
    with open(log_f_loc, 'r') as log_f:
        for line in log_f:
            done_list.append(line.strip())
    return done_list


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
