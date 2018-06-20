import subprocess

for i in range(19, 20):
    command = "tabix -h ../1-00801.hg38.trio.vcf.gz chr" + i + " > 1-00801_chr" + i + ".vcf";
    subprocess.call(command, shell=True);
