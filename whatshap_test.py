import subprocess
from WhatshapData import WhatshapData

patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

for ID in patientID:
    pedFileName = "hpc/users/seidea02/longreadclustersequencing/data/" + ID + ".ped";
    child_input = "hpc/users/seidea02/www/PacbioProject/PacbioHg38Bams/" + ID + "_edit_sorted.bam ";
    dad_input = "hpc/users/seidea02/www/PacbioProject/IlluminaHg38Bams/" + ID + "-02.hg38.dedup.clean.recal.cram ";
    mom_input = "hpc/users/seidea02/www/PacbioProject/IlluminaHg38Bams/" + ID + "-01.hg38.dedup.clean.recal.cram ";
    for i in range(1, 23):
        phaseFileName = ID + "_chr" + str(i) + "_phased.vcf ";
        inputVCF = "hpc/users/seidea02/www/PacbioProject/DNV_calls/VCF/TrioVCF/" + ID + "/" + ID + "_chr" + str(i) + ".vcf.gz ";
        whatshap_object = WhatshapData(ID, pedFileName, phaseFileName, inputVCF, child_input, dad_input, mom_input);
        print("========WhatshapData object created for patient family " + ID + " for chromosome " + str(i));
        command = whatshap_object.cmd();
        print("--------Running whatshap for patient family " + ID + " on chromosome " + str(i));

        print(command);
        #subprocess.call(command, shell=True);
