import subprocess as sp
import multiprocessing as mp

# patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

def whatshap(num):
    vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
    vcf_filename += "DNV_calls/VCF/TrioVCF/1-05846/1-05846_chr" + str(num) + ".vcf.gz";
    bam_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/PacbioHg38Bams/";
    bam_filename += "1-05846_edit_sorted.bam";
    command = "whatshap phase --sample=1-05846 --ignore-read-groups";
    command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
    command += " -o 1-05846_no_indels/1-05846_chr" + str(num) + "_phased.vcf ";
    command += vcf_filename + " " + bam_filename;
    sp.call(command, shell=True);
    print("======Sucessfully ran whatshap for 1-05846 on chromosome " + str(num));

def get_gtf(num):
    command = "whatshap stats --gtf=1-05846_chr" + str(num) + "phased.gtf 1-05846_no_indels/1-05846_chr" + str(num) + "_phased.vcf";
    sp.call(command, shell=True);
    print("------Succesfully created gtf for 1-05846 on chromosome " + str(num));

if __name__ == '__main__':
  pool = mp.Pool(processes=5);
  pool.map(whatshap, range(1,23));
  pool = mp.Pool(processes=5);
  pool.map(get_gtf, range(1,23));
