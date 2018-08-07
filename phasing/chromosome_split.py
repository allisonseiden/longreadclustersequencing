import subprocess
import multiprocessing as mp

# Script to split VCF files into 22 separate VCF files by chromosome number

# patientID = ["1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

# for ID in patientID:
#     cd = "cd /hpc/users/seidea02/www/PacbioProject/DNV_calls/VCF/TrioVCF/" + ID;
#     subprocess.call(cd, shell=True);
#     for i in range(1, 23):
#         filename = ID + "_chr" + str(i) + ".vcf";
#         split = "tabix -h /hpc/users/seidea02/www/PacbioProject/DNV_calls/VCF/TrioVCF/" + ID + ".hg38.trio.vcf.gz chr" + str(i) + " > " + filename;
#         subprocess.call(split, shell=True);
#         compress = "bgzip " + filename;
#         subprocess.call(compress, shell=True);
#         index = "tabix -p vcf " + filename + ".gz";
#         subprocess.call(index, shell=True);

def split_compress_index(num):
    filename = "Illumina_WGS_chr" + str(i) + ".vcf";
    split = "tabix -h /sc/orga/projects/chdiTrios/GMKF_WGS_Trios_Dec_2017/GMKF_Seidman_CongenitalHeartDisease_WGS.vcf.gz chr" + str(i) + " > " + filename;
    subprocess.call(split, shell=True);
    compress = "bgzip " + filename;
    subprocess.call(compress, shell=True);
    index = "tabix -p vcf " + filename + ".gz";
    subprocess.call(index, shell=True);

if __name__ == '__main__':
    pool = mp.Pool(processes=5);
    pool.map(split_compress_index, range(2, 23))
