import subprocess as sp
import multiprocessing as mp

patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05673", "1-05846"];

def get_gtf(ID):
    for i in range(1,23):
        command = "whatshap stats --gtf=" + ID + "_chr" + str(i) + "phased.gtf " + ID;
        command += "/" + ID + "_chr" + str(i) + "_phased.vcf";
        sp.call(command, shell=True);
        print("Created gtf for " + ID + " chromosome " + str(i));


if __name__ == '__main__':
    pool = mp.Pool(processes=5);
    pool.map(get_gtf, patientID)
