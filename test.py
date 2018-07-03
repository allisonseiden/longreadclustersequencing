from PhasedData import PhasedData


patient = PhasedData("1-00801");
patient.create_vcf_dictionary();
patient.create_dnvs_dictionary();
patient.fill_bounds_dictionary();
patient.find_variants_for_phasing();
patient.assign_to_parent();
# patient.convert_to_dataframe();
print(patient.phased_to_parent);
