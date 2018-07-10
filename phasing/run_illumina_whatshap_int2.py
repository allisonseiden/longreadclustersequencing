import subprocess as sp
import multiprocessing as mp

patientID = ["1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

def illumina_whatshap(ID):
    for i in range(1, 23):
        vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
        vcf_filename += "DNV_calls/VCF/TrioVCF/" + ID + "/" + ID + "_chr";
        vcf_filename += str(i) + ".vcf.gz"
        bam_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
        bam_filename += "IlluminaHg38Bams/" + ID + ".hg38.dedip.clean.recal.cram";
        command = "whatshap phase --sample=" + ID + " --ignore-read-groups";
        command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
        command += " --indels -o " + ID + "_illumina/" + ID + "_chr" + str(i) + "_phased.vcf ";
        command += vcf_filename + " " + bam_filename;
        sp.call(command, shell=True);
        print("======Sucessfully ran whatshap for " + ID + " on chromosome " + str(i));

if __name__ == '__main__':
  pool = mp.Pool(processes=5);
  pool.map(illumina_whatshap, patientID); 
