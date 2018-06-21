import subprocess
from WhatshapData import WhatshapData
import multiprocessing as mp

patientID = ["1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

#for ID in patientID:
def whatshap(ID):
    pedFileName = "/hpc/users/seidea02/longreadclustersequencing/data/" + ID + ".ped";
    child_input = "/hpc/users/seidea02/www/PacbioProject/PacbioHg38Bams/" + ID + "_edit_sorted.bam ";
    dad_input = "/hpc/users/seidea02/www/PacbioProject/IlluminaHg38Bams/" + ID + "-02.hg38.dedup.clean.recal.cram ";
    mom_input = "/hpc/users/seidea02/www/PacbioProject/IlluminaHg38Bams/" + ID + "-01.hg38.dedup.clean.recal.cram ";
    for i in range(1, 23):
        chrom = "chr{0}".format(i);
        phaseFileName = ID + "_" + chrom + "_phased.vcf ";
        inputVCF = "/hpc/users/seidea02/www/PacbioProject/DNV_calls/VCF/TrioVCF/" + ID + "/" + ID + "_" + chrom + ".vcf.gz ";
        whatshap_object = WhatshapData(ID, pedFileName, phaseFileName, inputVCF, child_input, dad_input, mom_input);
        print("========WhatshapData object created for patient family " + ID + " on chromosome " + chrom);
        command = whatshap_object.cmd();
        print("--------Running whatshap for patient family " + ID + " on chromosome " + chrom);
        subprocess.call(command, shell=True);
        print("********Succesfully ran whatshap for patient family " + ID + " on chromosome " + chrom);

if __name__ == '__main__':
    pool = mp.Pool(processes=6);
    pool.map(whatshap, patientID);
