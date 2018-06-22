import subprocess as sp

for i in range(1,15):
    vcf_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
    vcf_filename += "DNV_calls/VCF/TrioVCF/1-01019/1-01019_chr" + str(i) + ".vcf.gz";
    bam_filename = "/sc/orga/projects/chdiTrios/WGS_Combined_2017/PacbioProject/";
    bam_filename += "PacbioHg38Bams/1-01019_edit_sorted.bam";
    command = "whatshap phase --sample=1-01019 --ignore-read-groups";
    command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
    command += " --indels -o 1-01019/1-01019_chr" + str(i) + "_phased.vcf ";
    command += vcf_filename + " " + bam_filename;
    sp.call(command, shell=True);
    print("======Sucessfully ran whatshap for 1-01019 on chromosome " + str(i));
