#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions called in multiple files

:Authors: Allison Seiden, Felix Richter
:Date: 2018-08-19
:Copyright: 2018, Allison Seiden
:License: CC BY-SA

"""

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
