from PhasedData import PhasedData

# patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460",
#                 "1-04537", "1-05443", "1-05673", "1-05846"];

patientIDs = ["1-00801"];
phased_data_objects = [];

for ID in patientIDs:
    phased_data_objects.append(PhasedData(ID));

for patient_object in phased_data_objects:
    patient_object.create_vcf_dictionary();
    patient_object.create_dnvs_dictionary();
    patient_object.fill_bounds_dictionary();
    patient_object.find_variants_for_phasing();
    patient_object.assign_to_parent();
    patient_object.count_mom_and_dad();
    #patient_object.print_mom_and_dad_count();
    patient_object.convert_to_dataframe();
    #print(patient_object.num_each_parent);
