# Allison Seiden
# PacBio Project, creating class for Whatshap data to run shell commands

class WhatshapData:
    def __init__(self, patientID, pedFileName, phaseFileName, inputVCF, child_input, dad_input, mom_input):
        self.patient = patientID;
        self.ped = pedFileName;
        self.phasedName = phaseFileName;
        self.variantFileName = inputVCF;
        self.childseq = child_input;
        self.dadseq = dad_input;
        self.momseq = mom_input;
        self.command = "";

    def cmd(self):
        self.command = "whatshap phase --ped " + self.ped;
        self.command += " --distrust-genotypes --indels";
        self.command += " --reference /sc/orga/projects/chdiTrios/Felix/dbs/hg38.fa";
        self.command += " --include-homozygous -o " + self.phasedName + self.variantFileName;
        self.command += self.childseq + self.dadseq + self.momseq;
        return self.command;
