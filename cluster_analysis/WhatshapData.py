# Allison Seiden
# PacBio Project, creating class for Whatshap data to run shell commands

class WhatshapData(object):
    def __init__(self, patientID, pedFileName, phaseFileName, inputVCF, inputBAM):
        self.patient = patientID;
        self.ped = pedFileName;
        self.phasedName = phaseFileName;
        self.variantFileName = inputVCF;
        self.sequenceFileName = inputBAM;
        self.command = "";

    def cmd(self):
        self.command = "whatshap phase --ped " + self. ped + "--ignore-read-groups";
        self.command += " --distrust-genotypes --indels";
        self.command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
        self.command += " --include-homozygous -o " + self.phasedName;
        self.command += self.variantFileName + self.sequenceFileName;
        return self.command;
