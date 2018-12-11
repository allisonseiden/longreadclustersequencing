"""Meta-analysis of p-values

:Authors: Allison Seiden, Felix Richter
:Date: 2018-12-06
:Copyright: 2018, Allison Seiden, Felix Richter
:License: CC BY-SA

cd /Users/felixrichter/Dropbox/PhD/longreadclustersequencing/
python3

Sources
P-value meta-analysis: https://en.wikipedia.org/wiki/Fisher%27s_method
chi-squared to p-value: https://stackoverflow.com/a/11728072/10688049

"""

from scipy.stats import chi2
import numpy as np


# obtain pacbio_p_dict and ilmn_p_dict by copy-pasting from indels_illumina.py
# and create_indel_graphs.py

indel_list = ['HR', 'CCC', 'non-CCC', 'All']

pacbio_p_dict = {
    'Dad': {'HR': 0.7133596926454755, 'CCC': 0.008479624439501354,
            'non-CCC': 0.6950407935173517, 'All': 0.01623275487901527},
    'Mom': {'HR': 0.8673678868409915, 'CCC': 0.9992962802747138,
            'non-CCC': 0.1672429171043177, 'All': 0.34310231278041914}}

ilmn_p_dict = {
    'Dad': {'HR': 0.054525301895685746, 'CCC': 0.009148155194113357,
            'non-CCC': 0.055790766080513504, 'All': 0.0006536217962735643},
    'Mom': {'HR': 0.12522561889864206, 'CCC': 0.4605803210807925,
            'non-CCC': 0.3236462351108671, 'All': 0.1350704184236487}}

# 2k = degrees of freedom where k is number of tests
df = 2*2
for parent in ['Dad', 'Mom']:
    print(parent)
    for i_class in indel_list:
        print(i_class)
        pb_p = -2*np.log(pacbio_p_dict[parent][i_class])
        ilmn_p = -2*np.log(ilmn_p_dict[parent][i_class])
        x_i = pb_p + ilmn_p
        print(1 - chi2.cdf(x_i, df))
