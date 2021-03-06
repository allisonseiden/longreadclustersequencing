import subprocess

""" Script to sort and index edited BAM files once ascii characters are removed/
    read groups are added """

patientID = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460", "1-04537", "1-05443", "1-05673", "1-05846"];

for ID in patientID:
    sort = "samtools sort -o " + ID + "_edit_sorted.bam -O bam " + ID + "_edit.sam";
    index = "samtools index " + ID + "_edit_sorted.bam";
    subprocess.call(sort, shell=True);
    print("============" + ID + " sorted");
    subprocess.call(index, shell=True);
    print("============" + ID + " indexed");
