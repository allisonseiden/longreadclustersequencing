from PhasedData import PhasedData

patient_100801_indels = PhasedData("1-00801");
patient_100801_indels.create_vcf_dictionary();
patient_100801_indels.create_dnvs_dictionary();
patient_100801_indels.fill_bounds_dictionary();
patient_100801_indels.find_variants_for_phasing();
patient_100801_indels.assign_to_parent();


patient_100801_noindels = PhasedData("1-00801");
patient_100801_noindels.create_vcf_dictionary();
patient_100801_noindels.create_dnvs_dictionary();
patient_100801_noindels.fill_bounds_dictionary();
patient_100801_noindels.find_variants_for_phasing();
patient_100801_noindels.assign_to_parent();

for i in range(1,23):
    print('======With indels for chromosome chr' + str(i));
    print(patient_100801_indels.phased_to_parent["chr{0}".format(i)]);
    print('\n');
    print('------Without indels for chromosome chr' + str(i));
    print(patient_100801_noindels.phased_to_parent["chr{0}".format(i)]);
    print('\n\n');
