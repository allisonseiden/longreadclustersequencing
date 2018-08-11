import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb


sb.set_context("paper");
sb.set_style("ticks", {'font.family' : ['serif']});
# dnv_data = pd.read_csv('dnv_parent_percentages.txt', sep='\t', engine='python');
# dnv_data.set_index(['Indel_Class'], inplace=True);
# indel_list = ['HR', 'CCC', 'non-CCC'];
# length = len(indel_list);
#
# f, axarr = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10, 3));
# y_pos = [0.2625, 0.7375];
# plt.yticks([]);
#
#
# for i in range(length):
#     curr_df = dnv_data.loc[indel_list[i]];
#     curr_df.reset_index(inplace=True);
#     fraction_list = curr_df['Fraction'].tolist();
#     if i == 1:
#         dad = axarr[i].barh(y_pos[0], fraction_list[0], color="#352846", left=0.005, height=0.425, label='Father');
#         mom = axarr[i].barh(y_pos[1], fraction_list[1], color="#E2AE38", left=0.005, height=0.425, label='Mother');
#         axarr[i].set_title(indel_list[i]);
#     else:
#         axarr[i].barh(y_pos, fraction_list, color=["#352846",'#E2AE38'], left=[0.005, 0.005], height=0.425);
#         axarr[i].set_title(indel_list[i]);
#
#
# f.text(0.5, 0.06, 'Fraction', ha='center', va='center');
# f.text(0.02, 0.5, 'Indel Class', ha='center', va='center', rotation='vertical');
# f.tight_layout(pad=2.5, w_pad=1.0, h_pad=1.0);
# f.legend(handles=[dad, mom],loc=(0.776, 0.02), ncol=2, borderaxespad=0.);
#
# plt.savefig("mom_dad_comparison_bar_chart.png");

# -------------------- end of mom and dad bar mom_dad_comparison_bar_chart

# patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];
# indel_list = ['HR', 'CCC', 'non-CCC'];
# length = len(indel_list);
#
# parental_age = pd.read_csv('/Users/allisonseiden/Documents/longreadclustersequencing/data/parental_age_at_conception.txt', sep='\t', engine='python');
# parental_age.set_index(['ID'], inplace=True);
# parental_totals = pd.read_csv('/Users/allisonseiden/Documents/longreadclustersequencing/indel_analysis/dnv_parent_totals_by_ID.txt', sep='\t', engine='python');
# parental_totals.set_index(['ID'], inplace=True);
# dnv_data = parental_age.join(parental_totals, how='left');
#
#
# paternal_data = dnv_data[dnv_data.Parent == 'Father'];
# maternal_data = dnv_data[dnv_data.Parent == 'Mother'];
#
# f, axarr = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(10, 5));
# plt.xlim((0, 80));
# plt.ylim((-5, 10));
# for i in range(length):
#     # paternal_mut_df = paternal_data[paternal_data.Indel_Class == indel_list[i]];
#     # maternal_mut_df = maternal_data[maternal_data.Indel_Class == indel_list[i]];
#     age_list_dad = paternal_data['Paternal_age_at_conception'].tolist();
#     age_list_mom = maternal_data['Maternal_age_at_conception'].tolist();
#     total_list_dad = paternal_data[indel_list[i]].tolist();
#     total_list_mom = maternal_data[indel_list[i]].tolist();
#     dad = sb.regplot(x=age_list_dad, y=total_list_dad, scatter=True, color="#352846", ax=axarr[i], label='Father');
#     mom = sb.regplot(x=age_list_mom, y=total_list_mom, scatter=True, color="#E2AE38", ax=axarr[i], label='Mother');
#     axarr[i].set_title(indel_list[i]);
#
#
# f.text(0.02, 0.5, 'Total DNVs', ha='center', va='center', rotation='vertical');
# f.text(0.5, 0.035, 'Age of Parent at Conception (yr)', ha='center');
# f.tight_layout(pad=2.5, w_pad=1.0, h_pad=1.0);
# f.legend(labels=['Father', 'Mother'], ncol=2, loc=(0.772, 0.017));
# plt.savefig('/Users/allisonseiden/Documents/longreadclustersequencing/graphs/mom_dad_age_indel_totals_scatter.png');

#---------------- end of age scatter plot

indel_df = pd.read_table('/Users/allisonseiden/Documents/longreadclustersequencing/indel_analysis/all_indel_info.txt', sep='\t', engine='python');
length = indel_df.shape[0];

copy_count_df = indel_df[indel_df['CCC'] == 1];

print(copy_count_df);

# print(indel_df);
