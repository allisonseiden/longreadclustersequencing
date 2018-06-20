import subprocess

for i in range(1, 19):
    split = "tabix -h ../1-00801.hg38.trio.vcf.gz chr" + str(i) + " > 1-00801_chr" + str(i) + ".vcf";
    subprocess.call(split, shell=True);

    compress = "bgzip 1-00801_chr" + str(i) + ".vcf";
    subprocess.call(compress, shell=True);

    index = "tabix -p vcf 1-00801_chr" + str(i) + ".vcf.gz";
    subprocess.call(index, shell=True);
