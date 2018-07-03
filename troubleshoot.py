from PhasedData import PhasedData

patient_2 = PhasedData('1-01019');
patient_2.create_vcf_dictionary();
patient_2.create_dnvs_dictionary();
patient_2.fill_bounds_dictionary();
patient_2.find_variants_for_phasing();
patient_2.assign_to_parent();
patient_2.convert_to_dataframe();
print(patient_2.parent_df);
