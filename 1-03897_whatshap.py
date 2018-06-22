import subprocess as sp

for i in range(2,23):
    vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
    vcf_filename += "DNV_calls/VCF/TrioVCF/1-03897/1-03897_chr" + str(i) + ".vcf.gz";
    bam_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
    bam_filename += "PacbioHg38Bams/1-03897_edit_sorted.bam";
    command = "whatshap phase --sample=1-03897 --ignore-read-groups";
    command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
    command += " --indels -o 1-03897/1-03897_chr" + str(i) + "_phased.vcf ";
    command += vcf_filename + " " + bam_filename;
    sp.call(command, shell=True);
    print("======Sucessfully ran whatshap for 1-03897 on chromosome " + str(i));
