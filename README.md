# DNV phasing with long read technology

- Authors: Allison Seiden, Felix Richter


### Project overview

1. Use 5x PacBio whole genome sequencing (WGS) to phase de novo variants (DNVs) previously identified with 30x Illumina WGS.
    - Python 3 wrappers around WhatsHap to phase variants using PacBio (no indels) and/or Illumina WGS
    - Python 3 programs to assign DNVs to parent-of-origin based on informative variants (i.e., those uniquely inherited from mom or dad) on the same reads
2. Analyze properties of these DNVs, overall and between the following categories
    - Phased vs unphased
    - Illumina vs Illumina+PacBio phasing
    - Maternal vs paternal
    - Within clusters
    - As a function of parental age



### Dependencies
- Python 3
- [WhatsHap](whatshap.readthedocs.io)
- [tabix](www.htslib.org/doc/tabix.html) and [samtools](http://www.htslib.org/doc/samtools.html)


### Directory descriptions

- **cluster_analysis:** code for analyzing DNV clusters
- **data:** pedigree files for all trios
- **IGV:** salient IGV plots
- **notes:** lab notebook of day-to-day commands
- **phasing:** core programs to run WhatsHap and assign variants to parent-of-origin or classify as unphased
- **phasing_analysis:** scripts and programs for downstream DNV analyses



