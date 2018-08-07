""" Script to run Whatshap on all patient IDs without indel flag using multiprocessing
"""

import subprocess as sp
import multiprocessing as mp

patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

def whatshap(ID):
    for i in range(1, 23):
        vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
        vcf_filename += "DNV_calls/VCF/TrioVCF/" + ID + "/" + ID "_chr" + str(num) + ".vcf.gz";
        bam_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/PacbioHg38Bams/";
        bam_filename += ID + "_edit_sorted.bam";
        command = "whatshap phase --sample=" + ID + " --ignore-read-groups";
        command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
        command += " -o " + ID + "_no_indels/" + ID + "_chr" + str(num) + "_phased.vcf ";
        command += vcf_filename + " " + bam_filename;
        sp.call(command, shell=True);
        print("======Sucessfully ran whatshap for" + ID + " on chromosome " + str(num));

if __name__ == '__main__':
  pool = mp.Pool(processes=5);
  pool.map(whatshap, patientID);
