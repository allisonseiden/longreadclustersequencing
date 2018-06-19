import subprocess

subprocess.call("samtools sort -o 1-04460_edit_sorted.bam -O bam 1-04460_edit.bam", shell=True);
subprocess.call("samtools index 1-04460_edit_sorted.bam", shell=True);

subprocess.call("samtools sort -o 1-04537_edit_sorted.bam -O bam 1-04537_edit.bam", shell=True);
subprocess.call("samtools index 1-04537_edit_sorted.bam", shell=True);

subprocess.call("samtools sort -o 1-05443_edit_sorted.bam -O bam 1-05443_edit.bam", shell=True);
subprocess.call("samtools index 1-05443_edit_sorted.bam", shell=True);

subprocess.call("samtools sort -o 1-05673_edit_sorted.bam -O bam 1-05673_edit.bam", shell=True);
subprocess.call("samtools index 1-05673_edit_sorted.bam", shell=True);

subprocess.call("samtools sort -o 1-05846_edit_sorted.bam -O bam 1-05846_edit.bam", shell=True);
subprocess.call("samtools index 1-05846_edit_sorted.bam", shell=True);
