# **Annotating Metagenomically Assembled Bacteriophage from a Unique Ecological System using Protein Structure Prediction and Structure Homology Search**

Data was collected using Oxford Nanopore MinION Platform (Primarily LSK-110, R10.3/R9.4.1 Flowcells)
Raw reads were uploaded to the ENA [(PRJEB49151)](https://www.ebi.ac.uk/ena/browser/view/PRJEB49151)

**10 Samples have been included in the analysis so far:**

- GAC0 (2019)
- GAC1 (Jan 27, 20)
- GAC2 (Feb 14, 20)
- GAC3 (March 2, 20)
- GAC4 (Oct6, 20)
- GAC5(May 3, 21)
- GAC6 (May 7, 21)
- GAC7 (May 14, 21)
- GAC8 (Aug 10, 21)
- GAC9 (Apr 21, 22)

## Methods/Workflows

### [Assembly, Polishing and QC](./POLISHING-QC.md)

Contains steps used to assemble circular contigs and polish them. The polishing and coverage steps are also included as a [snakemake workflow](./snakemake_workflow).

### [Plasmid and Phage Analysis](./PLASMID_PHAGE_ANALYSIS.md)

Outlines methods to analyze and collect data on small circular contigs from this dataset.
