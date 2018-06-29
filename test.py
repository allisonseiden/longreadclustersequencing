from PhasedData import PhasedData


patient_100801 = PhasedData("1-00801");
patient_100801.create_vcf_dictionary();
patient_100801.create_dnvs_dictionary();
patient_100801.fill_bounds_dictionary();
patient_100801.find_variants_for_phasing_chr('chr22');
print(patient_100801.to_phase['chr22']);
patient_100801.assign_to_parent('chr22');
print(patient_100801.phased_to_parent);
