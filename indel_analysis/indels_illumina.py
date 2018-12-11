"""Analyze indels.

:Authors: Allison Seiden, Felix Richter
:Date: 2018-12-06
:Copyright: 2018, Allison Seiden, Felix Richter
:License: CC BY-SA

cd /Users/felixrichter/Dropbox/PhD/longreadclustersequencing/
python3

"""

import glob
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sb

# home_dir = '/hpc/users/richtf01/longreadclustersequencing/'
home_dir = '/Users/felixrichter/Dropbox/PhD/longreadclustersequencing/'

# load big phasing dataframe
dnv_f = home_dir + 'phasing_analysis/phasing_analysis_df_ilmn_2018_12_06.txt'
dnv_df = pd.read_table(dnv_f, sep='\t')
dnv_df.rename(columns={'Location': 'End'}, inplace=True)
dnv_df.head()
dnv_df.shape
dnv_df.Mom.isnull().sum()


# load indel info
indel_iter = glob.iglob(home_dir + 'data/gmkf2/sorting_hat_out/1-*.txt')
indel_list = []
for indel_f in indel_iter:
    indel_df = pd.read_table(indel_f, sep='\t')
    indel_list.append(indel_df)


indel_df = pd.concat(indel_list, ignore_index=True)
indel_df.head()
indel_df.shape

# 'Start', #
index_cols = ['ID', 'Chrom', 'End', 'Ref', 'Alt']
dnv_df.set_index(index_cols, inplace=True)
indel_df.set_index(index_cols, inplace=True)
indel_full_df = dnv_df.join(indel_df, how='inner')
indel_full_df.head()
# should have 2136 rows
indel_full_df.shape
indel_full_df.Mom.isnull().sum()

out_f = home_dir + 'phasing_analysis/indels_df_ilmn_2018_12_08.txt'
# indel_full_df.to_csv(out_f, sep='\t')


"""Plot parental age vs absolute indel count.

import pandas as pd
import numpy as np

"""

non_na_ids = dnv_df.loc[~dnv_df.Mom.isnull()].index.get_level_values(
    'ID').unique().tolist()
# indel_full_df = indel_full_df.loc[~indel_full_df.Mom.isnull()]
pacbio_ids = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389',
              '1-04460', '1-04537', '1-05443', '1-05673', '1-05846']

# patientIDs: get list of all patient IDs (including those w/o indels)
# but only those that were processed

indel_list = ['HR', 'CCC', 'non-CCC']
# indel_len = len(indel_list)
age_loc = home_dir + 'data/parental_age_at_conception_allPCGC_nonNA.txt'

age_df = pd.read_csv(age_loc, sep='\t', engine='python')
age_df.set_index(['ID'], inplace=True)
age_df = age_df.loc[non_na_ids]
# check if any of the PacBio IDs are among the GMKF2 IDs:
[i for i in pacbio_ids if i in age_df.index.values.tolist()]
# should be smae length/number of rows:
len(non_na_ids)
age_df.shape

# select for NAs, otherwise phased_indels includes NA indels
indel_nonNa_df = indel_full_df.loc[~indel_full_df.Mom.isnull()]
phased_indels = (indel_nonNa_df.Mom != 0) | (indel_nonNa_df.Dad != 0)
indel_phased_df = indel_nonNa_df.loc[phased_indels]

# overall stats on phasing
indel_nonNa_df.shape  # 1535
indel_phased_df.shape  # 306
indel_phased_df.Mom.sum()  # 82
indel_phased_df.Dad.sum()  # 224


def get_indel_ct_dict(id_i, id_df, indel_list):
    """Get counts per indel type for the subgroup."""
    indel_ct_dict = dict(zip(indel_list, [0, 0, 0]))
    indel_ct_dict['ID'] = id_i
    # print(id_df)
    # print(id_df['Indel_Class'])
    for indel_i in indel_list:
        # print(id_df['Indel_Class'].str.match(indel_i).sum())
        indel_ct = id_df['Indel_Class'].str.match(indel_i).sum()
        indel_ct_dict[indel_i] = indel_ct
    return indel_ct_dict


def get_indel_df(indel_phased_df, indel_list, age_df):
    """Get the dataframe of indel counts per ID."""
    id_group = indel_phased_df.groupby(
        [indel_phased_df.index.get_level_values('ID')])
    # get per-ID indels per class
    # count = 0
    indel_ct_list = []
    for id_i, id_df in id_group:
        # print(id_i)
        indel_ct_dict = get_indel_ct_dict(id_i, id_df, indel_list)
        # print(pd.crosstab(id_df['Indel_Class'], columns='count'))
        # print(indel_ct_dict)
        indel_ct_list.append(indel_ct_dict)
        # count += 1
        # if count > 10:
        #     break
    indel_class_df = pd.DataFrame(indel_ct_list)
    indel_class_df.set_index(['ID'], inplace=True)
    indel_class_df = age_df.join(indel_class_df, how='left').fillna(0)
    return indel_class_df


indel_df_dad = indel_phased_df.loc[indel_phased_df.Dad == 1]
indel_df_dad = get_indel_df(indel_df_dad, indel_list, age_df)

indel_df_mom = indel_phased_df.loc[indel_phased_df.Mom == 1]
indel_df_mom = get_indel_df(indel_df_mom, indel_list, age_df)

for indel_i in indel_list:
    print('{} dad ct: {}, mom ct: {}'.format(
        indel_i,  indel_df_dad[indel_i].sum(), indel_df_mom[indel_i].sum()))

"""Actually plot."""

sb.set_context("paper")
sb.set_style("ticks", {'font.family': ['serif']})

f, axarr = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10, 5))
plt.xlim((0, 80))
plt.ylim((-0.25, 3.25))
plt.yticks([i for i in range(4)])
for i_class in indel_list:
    print(i_class)
    i = indel_list.index(i_class)
    age_list_dad = indel_df_dad['Paternal_age_at_conception']  # .tolist()
    age_list_mom = indel_df_mom['Maternal_age_at_conception']
    total_list_dad = indel_df_dad[i_class]  # .tolist()
    total_list_mom = indel_df_mom[i_class]
    # print(age_list_dad, total_list_dad)
    dad = sb.regplot(x=age_list_dad, y=total_list_dad, scatter=True,
                     color="#352846", ax=axarr[i], label='Father')
    mom = sb.regplot(x=age_list_mom, y=total_list_mom, scatter=True,
                     color="#E2AE38", ax=axarr[i], label='Mother')
    axarr[i].set_title(i_class)
    axarr[i].xaxis.label.set_visible(False)


f.text(0.02, 0.5, 'Total DNVs', ha='center', va='center', rotation='vertical')
f.text(0.5, 0.035, 'Age of Parent at Conception (yr)', ha='center')
f.tight_layout(pad=2.5, w_pad=1.0, h_pad=1.0)
# f.legend(labels=['Father', 'Mother'], ncol=2, loc=(0.772, 0.017))
out_f = home_dir + 'graphs/mom_dad_age_indel_totals_scatter_ilmn_2018_12.png'
plt.savefig(out_f)

indel_df_dad['All'] = (
    indel_df_dad['HR'] + indel_df_dad['CCC'] + indel_df_dad['non-CCC'])
indel_df_mom['All'] = (
    indel_df_mom['HR'] + indel_df_mom['CCC'] + indel_df_mom['non-CCC'])

indel_list.append('All')
ilmn_p_dict = {'Dad': {}, 'Mom': {}}
for i_class in indel_list:
    print(i_class)
    age_list_dad = indel_df_dad['Paternal_age_at_conception']  # .tolist()
    age_list_mom = indel_df_mom['Maternal_age_at_conception']
    total_list_dad = indel_df_dad[i_class]  # .tolist()
    total_list_mom = indel_df_mom[i_class]
    print('Dad totals: {}, correlation:'.format(sum(total_list_dad)))
    pearsonr(age_list_dad, total_list_dad)
    ilmn_p_dict['Dad'][i_class] = pearsonr(age_list_dad, total_list_dad)[1]
    print('Mom totals: {}, correlation:'.format(sum(total_list_mom)))
    pearsonr(age_list_mom, total_list_mom)
    ilmn_p_dict['Mom'][i_class] = pearsonr(age_list_mom, total_list_mom)[1]


#
