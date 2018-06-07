# Allison Seiden
# PacBio Project, creating class for Whatshap data to run shell commands

class WhatshapData(object):
    def __init__(self, patientID, phaseFileName, inputVCF, inputBAM):
        self.patient = patientID;
        self.phasedName = phaseFileName;
        self.variantFileName = inputVCF;
        self.sequenceFileName = inputBAM;

    def cmd(self):
        self.command = "whatshap phase --ignore-read-groups --distrust-genotypes"
        self.command += " --indels --reference (insert ref file here)"
        self.command += " --include-homozygous -o ${output}" + self.phasedName;
        self.command += " ${output}" + self.variantFileName;
        self.command += " ${output}" + self.sequenceFileName;
        print(self.command);


test = WhatshapData(12345, "phaseFile", "VCFfile", "BAMfile");
test.cmd();
