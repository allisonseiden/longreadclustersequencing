from PhasedData import PhasedData

patient_2 = PhasedData('1-01019');
patient_2.create_vcf_dictionary();
patient_2.create_dnvs_dictionary();
patient_2.fill_bounds_dictionary();
print(patient_2.bounds);
