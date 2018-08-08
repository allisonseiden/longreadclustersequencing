import pandas as pd
import subprocess as sp
import multiprocessing as mp

""" Script 2/2 to run Whatshap with indels flag using Illumina data for second
    5 patient IDs, uses multiprocessing """

df = pd.read_table('/hpc/users/seidea02/longreadclustersequencing/data/Illumina_crams_batch2.txt', names=['ID']);
patientID = df['ID'].tolist();

# patientID = ["1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

def illumina_whatshap(ID):
    mkdir = 'mkdir ' + ID;
    sp.call(mkdir, shell=True);
    cd = 'cd /sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/IlluminaWhatshapVCFs/Batch2/' + ID;
    sp.call(cd, shell=True);
    for i in range(1, 23):
        vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
        vcf_filename += "IlluminaVCFs/Illumina_WGS_chr" + str(i) + ".vcf.gz";
        bam_filename = "/sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_2017/CRAM/Batch2/";
        bam_filename += ID + ".cram";
        command = "whatshap phase --sample=" + ID + " --ignore-read-groups";
        command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
        command += " --indels -o " + ID + "/" + ID + "_chr" + str(i) + "_phased.vcf ";
        command += vcf_filename + " " + bam_filename;
        sp.call(command, shell=True);
        print("======Sucessfully ran whatshap for " + ID + " on chromosome " + str(i));
    sp.call('cd ..', shell=True);

if __name__ == '__main__':
  pool = mp.Pool(processes=5);
  pool.map(illumina_whatshap, patientID);
