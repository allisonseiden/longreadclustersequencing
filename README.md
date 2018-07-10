# DNV phasing with long read technology

- Authors: Allison Seiden, Felix Richter


### Project overview

1. Use 5x PacBio whole genome sequencing (WGS) to phase de novo variants (DNVs) previously identified with 30x Illumina WGS.
    - Phase PacBio only with WhatsHap (no indels)
    - Use custom Python 3 programs to assign DNVs to parent-of-origin based on informative variants (i.e., those uniquely inherited from mom or dad) on the same reads
2. Analyze properties of these DNVs, overall and between the following categories
    - Phased vs unphased
    - Illumina vs Illumina+PacBio phasing
    - Maternal vs paternal
    - Within clusters
    - As a function of parental age


