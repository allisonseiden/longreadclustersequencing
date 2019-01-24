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

### Core phasing pipeline

- **chromosome_split.py:** splits a VCF by chromosome in (so that Whatshap is run separately by chromosome)
- **whasthap_bsub.sh (or illumina_whatshap_int1.py):** Run whatshap with indels per chromosome
    - **pacbio_whatshap.py:** run phasing analysis for PacBio data
- **whasthap_output_check.py:** Check whatshap output file length and delete if not consistent with input (i.e., re-run)
- **clean_whatshap_vcf.py:** Remove rows without variants and move VCF to new, smaller directory
- **get_gtf.py:** Run the Whasthap GTF command to get phased haplotype blocks
- **get_ID_dataframes.py:** wrapper around the PhasedData.py class, creates a PhasedData object for every ID/chromosome
- Other files in the **phasing/** directory:
    - check_discontinuities.py: tests the discontinuities functions that are part of PhasedData.py
    - get_results.py: review the phasing results
    - sort_and_index.py: preprocess the PacBio BAM (which had weird ASCII characters that had to be removed)
    - split_trios.py: split the large single VCF by trio
    - utilsy.py: various utility functions used across scripts

### Analysis file descriptions

- PacBio (N=10) data results
    1. phasing_analysis/results_phasing/indels_df.txt: indels phased with heuristic/IGV approach
    2. indel_analysis/classified_indels.txt: sorting-hat output (contains indel classes and repeat track overlaps)
    3. indel_analysis/all_indel_info.txt: joins of (1) and (2)
- Illumina (N=308) data results
    1. phasing_analysis/results_phasing/indels_df_ilmn_2018_12_08.txt: contains *both* indel phasing and sorting hat output classifications
- counts_per_id_ilmn_pb_2019_01_22.txt: Joined counts per ID from Illumina and PacBio

