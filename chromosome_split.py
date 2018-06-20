import subprocess

patientID = ["1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

for ID in patientID:
    cd = "/hpc/users/seidea02/www/PacbioProject/DNV_calls/VCF/TrioVCF/" + ID;
    subprocess.call(cd, shell=True);
    for i in range(1, 23):
        filename = ID + "_chr" + str(i) + ".vcf";
        split = "tabix -h ../" + ID + ".hg38.trio.vcf.gz chr" + str(i) + " > " + filename;
        subprocess.call(split, shell=True);
        compress = "bgzip " + filename;
        subprocess.call(compress, shell=True);
        index = "tabix -p vcf " + filename + ".gz";
        subprocess.call(index, shell=True);
