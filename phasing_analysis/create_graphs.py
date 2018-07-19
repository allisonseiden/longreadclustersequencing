import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import wes


sb.set_context("paper");
sb.set_style("ticks", {"font.family" : ['calibri']})
# dnv_data = pd.read_csv('dnv_parent_percentages.txt', sep='\t', engine='python');
# dnv_data.set_index(['Mutational_Class'], inplace=True);
# mutation_list = ['C_A', 'C_T', 'C_G', 'T_A', 'T_C', 'T_G', 'CpG_TpG'];
# length = len(mutation_list);
#
# f, axarr = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(10, 3));
# y_pos = [0.2625, 0.7375];
# plt.yticks([]);
#
# for i in range(length):
#     curr_df = dnv_data.loc[mutation_list[i]];
#     curr_df.reset_index(inplace=True);
#     fraction_list = curr_df['Fraction'].tolist();
#     if i < 4:
#         axarr[0, i].barh(y_pos, fraction_list, color=["#352846", '#E2AE38'], left=[0.005, 0.005], height=0.425, align='center');
#         axarr[0, i].set_title(mutation_list[i][0] + '  >  ' + mutation_list[i][2]);
#     elif i == 4:
#         axarr[1, 0].barh(y_pos, fraction_list, color=["#352846", '#E2AE38'], left=[0.005, 0.005], height=0.425, align='center');
#         axarr[1, 0].set_title(mutation_list[i][0] + '  >  ' + mutation_list[i][2]);
#     elif i == 5:
#         axarr[1, 1].barh(y_pos, fraction_list, color=["#352846", '#E2AE38'], left=[0.005, 0.005], height=0.425, align='center');
#         axarr[1, 1].set_title(mutation_list[i][0] + '  >  ' + mutation_list[i][2]);
#     else:
#         dad = axarr[1, 2].barh(y_pos[0], fraction_list[0], color='#352846', left=0.005, height=0.425, label='Father');
#         mom = axarr[1, 2].barh(y_pos[1], fraction_list[1], color='#E2AE38', left=0.005, height=0.425, label='Mother');
#         axarr[1, 2].set_title('CpG  >  TpG');
# f.text(0.5, 0.06, 'Fraction', ha='center', va='center');
# f.text(0.02, 0.5, 'Mutational Class', ha='center', va='center', rotation='vertical');
# f.tight_layout(pad=2.5, w_pad=1.0, h_pad=1.0);
# f.legend(handles=[dad, mom],loc=(0.776, 0.02), ncol=2, borderaxespad=0.);
#
# plt.savefig("mom_dad_comparison_bar_chart.png");
#
# -------------------- end of mom and dad comparison bar chart
#
# patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];
# mutation_list = ['C_A', 'C_T', 'C_G', 'T_A', 'T_C', 'T_G', 'CpG_TpG'];
# length = len(mutation_list);
#
# parental_age = pd.read_csv('/Users/Seiden/Documents/Summer2018/longreadclustersequencing/data/parental_age_at_conception.txt', sep='\t', engine='python');
# parental_age.set_index(['ID'], inplace=True);
# parental_fractions = pd.read_csv('/Users/Seiden/Documents/Summer2018/longreadclustersequencing/phasing_analysis/dnv_parent_percentages_by_ID.txt', sep='\t', engine='python');
# parental_fractions.set_index(['ID'], inplace=True);
# dnv_data = parental_age.join(parental_fractions, how='left');
#
# paternal_data = dnv_data[dnv_data.Parent == 'Father'];
# maternal_data = dnv_data[dnv_data.Parent == 'Mother'];
#
# f, axarr = plt.subplots(2, 4, sharex=True, sharey=True, figsize=(10, 5));
# plt.xlim((0, 100));
# for i in range(length):
#     paternal_mut_df = paternal_data[paternal_data.Mutational_Class == mutation_list[i]];
#     maternal_mut_df = maternal_data[maternal_data.Mutational_Class == mutation_list[i]];
#     age_list_dad = paternal_mut_df['Paternal_age_at_conception'].tolist();
#     age_list_mom = maternal_mut_df['Maternal_age_at_conception'].tolist();
#     fraction_list_dad = paternal_mut_df['Fraction'].tolist();
#     fraction_list_mom = maternal_mut_df['Fraction'].tolist();
#     if i < 4:
#         dad = sb.regplot(x=age_list_dad, y=fraction_list_dad, scatter=True, color="#352846", ax=axarr[0, i], label='Father');
#         mom = sb.regplot(x=age_list_mom, y=fraction_list_mom, scatter=True, color="#E2AE38", ax=axarr[0, i], label='Mother');
#         axarr[0, i].set_title(mutation_list[i][0] + '  >  ' + mutation_list[i][2]);
#     elif i == 4:
#         dad = sb.regplot(x=age_list_dad, y=fraction_list_dad, scatter=True, color="#352846", ax=axarr[1, 0], label='Father');
#         mom = sb.regplot(x=age_list_mom, y=fraction_list_mom, scatter=True, color="#E2AE38", ax=axarr[1, 0], label='Mother');
#         axarr[1, 0].set_title(mutation_list[i][0] + '  >  ' + mutation_list[i][2]);
#     elif i == 5:
#         dad = sb.regplot(x=age_list_dad, y=fraction_list_dad, scatter=True, color="#352846", ax=axarr[1, 1], label='Father');
#         mom = sb.regplot(x=age_list_mom, y=fraction_list_mom, scatter=True, color="#E2AE38", ax=axarr[1, 1], label='Mother');
#         axarr[1, 1].set_title(mutation_list[i][0] + '  >  ' + mutation_list[i][2]);
#     else:
#         dad = sb.regplot(x=age_list_dad, y=fraction_list_dad, scatter=True, color="#352846", ax=axarr[1, 2], label='Father');
#         mom = sb.regplot(x=age_list_mom, y=fraction_list_mom, scatter=True, color="#E2AE38", ax=axarr[1, 2], label='Mother');
#         axarr[1, 2].set_title('CpG  >  TpG');
#
#
# f.text(0.02, 0.5, 'Fraction', ha='center', va='center', rotation='vertical');
# f.text(0.5, 0.035, 'Age of Parent at Conception (yr)', ha='center');
# f.tight_layout(pad=2.5, w_pad=1.0, h_pad=1.0);
# f.legend(labels=['Father', 'Mother'], ncol=2, loc=(0.772, 0.017));
# plt.savefig('/Users/Seiden/Documents/Summer2018/longreadclustersequencing/graphs/mom_dad_age_scatter.png');
#
# --------------------------- end of parental age scatter plot

pacbio_data = pd.read_table('all_pacbio_data.txt', sep='\t', engine='python');
# f, axarr = plt.subplots(2, 5, sharex=True, sharey=True);
# chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr7', 'chr8', 'chr9', 'chr10', 'chr12'];
# i = 0;
# for chr in chr_list:
#     chr_df = pacbio_data[pacbio_data.Chrom == chr];
#     chr_df_mom = chr_df[chr_df.Parent == 'Mother'];
#     loc_list = chr_df_mom['Location'].tolist();
#     if i < 5:
#         axarr[0, i].hist(loc_list, color="#E2AE38");
#         axarr[0, i].set_title(chr);
#     else:
#         axarr[1, i-5].hist(loc_list, color="#E2AE38");
#         axarr[1, i-5].set_title(chr);
#     i += 1;
# f.tight_layout();

f_2, axarr_2 = plt.subplots(3, 6, sharex=True, sharey=True);
for i in range(1, 19):
    chrom_name = 'chr' + str(i);
    chr_df = pacbio_data[pacbio_data.Chrom == chrom_name];
    chr_df_dad = chr_df[chr_df.Parent == 'Father'];
    loc_list = chr_df_dad['Location'].tolist();
    if i < 7:
        axarr_2[0, (i-1)].hist(loc_list, color="#352846");
        axarr_2[0, (i-1)].set_title(chrom_name);
    elif i < 13:
        axarr_2[1, i-7].hist(loc_list, color="#352846");
        axarr_2[1, i-7].set_title(chrom_name);
    else:
        axarr_2[2, i-13].hist(loc_list, color="#352846");
        axarr_2[2, i-13].set_title(chrom_name);
f_2.tight_layout();

plt.savefig('/Users/Seiden/Documents/Summer2018/longreadclustersequencing/graphs/chromosomal_regions_dad.png');
