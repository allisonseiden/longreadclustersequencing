import pandas as pd
import subprocess as sp
import multiprocessing as mp

""" Script 1/2 to run Whatshap with indels flag using Illumina data for first
    5 patient IDs, uses multiprocessing """

df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/Illumina_crams_batch1.txt', names=['ID']);
patientID = df['ID'].tolist();

# patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389"];

def illumina_whatshap(ID):
    mkdir = 'mkdir ' + ID;
    sp.call(mkdir, shell=True);
    cd = 'cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/IlluminaWhatshapVCFs/Batch1' + ID;
    sp.call(cd, shell=True);
    for i in range(1, 23):
        vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
        vcf_filename += "IlluminaVCFs/Illumina_WGS_chr" + str(i) + ".vcf.gz";
        bam_filename = "/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_2017/CRAM/Batch1/";
        bam_filename += ID + ".cram";
        command = "whatshap phase --sample=" + ID + " --ignore-read-groups";
        command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
        command += " --indels -o " + ID + "/" + ID + "_chr" + str(i) + "_phased.vcf ";
        command += vcf_filename + " " + bam_filename;
        sp.call(command, shell=True);
        print("======Sucessfully ran whatshap for " + ID + " on chromosome " + str(i));

if __name__ == '__main__':
  pool = mp.Pool(processes=5);
  pool.map(illumina_whatshap, patientID);
