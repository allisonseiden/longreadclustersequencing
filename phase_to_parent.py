from PhasedData import PhasedData;
import pandas as pd;
import numpy as np;

patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460",
                "1-04537", "1-05443", "1-05673", "1-05846"];

phased_data_objects = [];
complete_df = pd.DataFrame();

for ID in patientIDs:
    phased_data_objects.append(PhasedData(ID));

for patient_object in phased_data_objects:
    patient_object.create_vcf_dictionary();
    patient_object.create_dnvs_dictionary();
    patient_object.fill_bounds_dictionary();
    patient_object.find_variants_for_phasing();
    patient_object.assign_to_parent();
    patient_object.convert_to_dataframe();
    complete_df.append(patient_object.parent_df);
