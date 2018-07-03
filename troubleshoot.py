from PhasedData import PhasedData

# patient_2 = PhasedData('1-01019');
# patient_2.create_vcf_dictionary();
# patient_2.create_dnvs_dictionary();
# patient_2.fill_bounds_dictionary();
# patient_2.find_variants_for_phasing();
# patient_2.assign_to_parent();
# patient_2.convert_to_dataframe();
# print(patient_2.parent_df);

# patient_3 = PhasedData('1-03897');
# patient_3.create_vcf_dictionary();
# patient_3.create_dnvs_dictionary();
# patient_3.fill_bounds_dictionary();
# patient_3.find_variants_for_phasing();
# patient_3.assign_to_parent();
# patient_3.convert_to_dataframe();
# print(patient_3.parent_df);

# patient_4 = PhasedData('1-04190');
# patient_4.create_vcf_dictionary();
# patient_4.create_dnvs_dictionary();
# patient_4.fill_bounds_dictionary();
# patient_4.find_variants_for_phasing();
# patient_4.assign_to_parent();
# patient_4.convert_to_dataframe();
# print(patient_4.parent_df);


patient_5 = PhasedData('1-04389');
patient_5.create_vcf_dictionary();
patient_5.create_dnvs_dictionary();
patient_5.fill_bounds_dictionary();
patient_5.find_variants_for_phasing();
patient_5.assign_to_parent();
patient_5.convert_to_dataframe();
print(patient_5.parent_df);
