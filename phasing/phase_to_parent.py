from PhasedData import PhasedData;
import pandas as pd;
import numpy as np;
import multiprocessing as mp;

# removed "1-04389", "1-04537" from list because of errors
# patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04460", "1-05443",
#                 "1-05673", "1-05846"];

patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04460"];

phased_df_list = [];
phased_data_objects = [];

for ID in patientIDs:
    phased_data_objects.append(PhasedData(ID));
# # file = open('proband_phased_data.txt', 'w');
#
# # for ID in patientIDs:
# def run_phase_to_parent(ID):
#     patient = PhasedData(ID);
#     patient.create_vcf_dictionary();
#     patient.create_dnvs_dictionary();
#     patient.fill_bounds_dictionary();
#     patient.find_variants_for_phasing();
#     patient.assign_to_parent();
#     patient.convert_to_dataframe();
#     # complete_df.append(patient.parent_df);
#     # file.write(patient.parent_df);
#     phased_df.append(patient.parent_df);
#
# if __name__ == '__main__':
#     pool = mp.Pool(processes=5);
#     pool.map(run_phase_to_parent, patientIDs);
#
#     print(phased_df);
#     # file.close();
#     #bigdata = pd.concat(phased_df, ignore_index=True);
#     #print(bigdata);
#     #for patient in phased_data_objects:
#     #     # complete_df.append(patient.parent_df);
#     #     print(patient.parent_df);
#
#     #print(complete_df);



for patient_object in phased_data_objects:
    patient_object.create_vcf_dictionary();
    patient_object.create_dnvs_dictionary();
    patient_object.fill_bounds_dictionary();
    patient_object.find_variants_for_phasing(3);
    patient_object.assign_to_parent();
    patient_object.convert_to_dataframe();
    phased_df_list.append(patient_object.parent_df);

complete_df = pd.concat(phased_df_list, ignore_index=True);

complete_df = complete_df.groupby('ID').sum();
complete_df = complete_df.loc[:,['From Mom', 'From Dad', 'Troubleshoot', 'Unphased']];

complete_df.to_csv();