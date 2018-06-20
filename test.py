import subprocess

for i in range(19, 20):
    command = "tabix -h ../1-00801.hg38.trio.vcf.gz chr" + str(i) + " > 1-00801_chr" + str(i) + ".vcf";
    subprocess.call(command, shell=True);
