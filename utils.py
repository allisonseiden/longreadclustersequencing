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
