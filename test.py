from PhasedData import PhasedData


patient_100801 = PhasedData("1-00801");
patient_100801.create_vcf_dictionary();
patient_100801.create_dnvs_dictionary();
patient_100801.fill_bounds_dictionary();
patient_100801.find_variants_for_phasing();
patient_100801.assign_to_parent();
patient_100801.convert_to_dataframe();
print(patient_100801.parent_df);
