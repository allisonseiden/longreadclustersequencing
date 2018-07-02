from PhasedData import PhasedData;
import pandas as pd;
import numpy as np;
import multiprocessing as mp;

patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460",
                 "1-04537", "1-05443", "1-05673", "1-05846"];

phased_data_objects = [];
complete_df = pd.DataFrame();

# for ID in patientIDs:
def run_phase_to_parent(ID):
    patient = PhasedData(ID);
    patient.create_vcf_dictionary();
    patient.create_dnvs_dictionary();
    patient.fill_bounds_dictionary();
    patient.find_variants_for_phasing();
    patient.assign_to_parent();
    patient.convert_to_dataframe();
    phased_data_objects.append(patient);

if __name__ == '__main__':
    pool = mp.Pool(processes=3);
    pool.map(run_phase_to_parent, patientIDs);

    for patient in phased_data_objects:
        complete_df.append(patient.parent_df);

    print(complete_df);



# for patient_object in phased_data_objects:
#     patient_object.create_vcf_dictionary();
#     patient_object.create_dnvs_dictionary();
#     patient_object.fill_bounds_dictionary();
#     patient_object.find_variants_for_phasing();
#     patient_object.assign_to_parent();
#     patient_object.convert_to_dataframe();
#     complete_df.append(patient_object.parent_df);

# print(complete_df);
