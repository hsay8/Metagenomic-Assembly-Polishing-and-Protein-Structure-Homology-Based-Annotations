# Polishing and QC Steps

- [Polishing and QC Steps](#polishing-and-qc-steps)
  - [QC and Draft Assembly](#qc-and-draft-assembly)
    - [Basecalling (Guppy)](#basecalling-guppy)
    - [QC Plots of Basecalled Reads (pycoQC)](#qc-plots-of-basecalled-reads-pycoqc)
    - [Filtering Basecalled Reads (NanoFilt)](#filtering-basecalled-reads-nanofilt)
    - [Assembly (Flye)](#assembly-flye)
    - [Extract Assembled Circular Contigs (ACCs)](#extract-assembled-circular-contigs-accs)
  - [Polishing](#polishing)
    - [Minimap2](#minimap2)
    - [Gerenuq](#gerenuq)
    - [Round 1: Minipolish](#round-1-minipolish)
    - [Round 2: Medaka](#round-2-medaka)
    - [QC: Coverage and Completeness](#qc-coverage-and-completeness)
    - [Polished and Filtered Coverage](#polished-and-filtered-coverage)
    - [Completeness](#completeness)
    - [Running this on Secondary Assemblies](#running-this-on-secondary-assemblies)

---

Sequences were obtained using Nanopore's MinION platform
Raw reads on gru were compressed to tar.gz to be uploaded to the ENA **(Study ascession: PRJEB49151)**

Note: this is also available as a [Snakemake pipeline](snakemake_workflow/Snakefile).

**Directory Tree:**

```bash
#Set Working Directory
WKDIR=''
├── 2-nanofilt
│   └── basecalled_filtered.fastq
├── 3-assemblyhq
├── 4-minimap
├── 6-minipolish
├── 7-medaka
├── checkm
├── envs
├── final_assemblies
├── gcinfo
├── gtdbtk_taxonomy
├── logs
├── metaphlan
├── mosdepth
├── polished_depth
├── prokka
├── Snakefile
└── tools

#Read data
RAWREADS=./rawreads
BASECALLED=$WKDIR/1-basecalled

#Reads QC
QCPLOTS=$WKDIR/2-nanofilt/nanoplot

#Polishing
FILTREADS=$WKDIR/2-nanofilt/basecalled_filtered.falstq
ASSEMBLY=$WKDIR/3-assemblyhq
CIRC_CONTIGS=$WKDIR/3-assemblyhq/circularized
MAPPED=$WKDIR/4-minimap
GQ_FILT=$WKDIR/5-gerenuq
MINIPOL=$WKDIR/6-minipolish
MEDAKA=$WKDIR/7-medaka
FINAL_ASSEMBLY=$WKDIR/final_assemblies

#Coverage
DEPTH_POLISH=$WKDIR/mosdepth/polished
DEPTH_UNPOLISH=$WKDIR/mosdepth/unpolished
MAP_POLISHED=$WKDIR/polished_depth

#Assembly QC and Analysis
CHECKM=$WKDIR/checkm
GTDB=$WKDIR/gtdbtk_taxonomy
PROKKA=$WKDIR/prokka
ANVIO=$WKDIR/anvio

#Set path for scripts
GETCIRCULARIZED='scripts/extract-circularized-fasta.py'
GFATOOLS='scripts/gfatools/gfatools'
MINIASM_FORM='scripts/flyetominiasm.py'
```

Create directories (the rest will be created by Snakemake):

```bash
mkdir -p $WKDIR $WKDIR/0-rawreads $WKDIR/2-nanofilt $WKDIR/2-nanofilt/pycoQC_plots $WKDIR/3-assemblyhq $WKDIR/3-assemblyhq/circularized_assemblies $WKDIR/mosdepth $WKDIR/mosdepth/polished $WKDIR/mosdepth/unpolished $WKDIR/polished_mapped$WKDIR/prokka $WKDIR/prokka/16s_sequences $WKDIR/gcinfo $WKDIR/checkm $WKDIR/gtdbtk_taxonomy $WKDIR/final_assemblies
```

---

## QC and Draft Assembly

### Basecalling (Guppy)

Basecalled in super accuracy mode config, keep reads that have q-score > 7

```bash
nohup guppy_basecaller -i $RAWREADS -r -s $BASECALLED -c dna_r9.4.1_450bps_sup.cfg -x cuda:0 --min_qscore 7 &> $BASECALLED/basecalling.log &
```

### QC Plots of Basecalled Reads (pycoQC)

```bash
pycoQC -f $BASECALLED/sequencing_summary.txt -o $QCPLOTS/pycoQC_output.html
```

### Filtering Basecalled Reads (NanoFilt)

```bash
LENGTH=500
Q_SCORE=7

#Filter basecalls by length and quality into a new fastq
NanoFilt -l $LENGTH -q $Q_SCORE --headcrop 50 < $BASECALLED/basecalled_concat.fastq > $FILTREADS &
```

### Assembly (Flye)

When using super accuracy mode basecalls (error <5%), if the error rate distribution is mostly below threshold (viewed through the flye.log file, under overlap based coverage), nano-hq is fine. Otherwise use raw.

```bash
#first attempt with hq mode. this mode is optimized for error<5% (sup-reads)
nohup flye --nano-hq $FILTREADS -o $ASSEMBLY -t 112 --meta &> $ASSEMBLY/assembly.log &
```

Remove circularized contigs that are flagged by Flye as repetitive, and have estimated coverages less than 10.

```bash
#get list of circular contigs that are repetitive and have coverages of less than 10
awk '$3 < 10 && $4 == "Y" && $5 == "Y" {print $1}' $ASSEMBLY/assembly_info.txt > $ASSEMBLY/contigs_to_remove.txt

#remove them
while read -r contig; do rm $ASSEMBLY/circularized_assemblies/$contig; done < $ASSEMBLY/contigs_to_remove.txt
```

### Extract Assembled Circular Contigs (ACCs)

Extract the fasta of ACCs above a specified size threshold (bases) into the $ASSEMBLY/circularized_assemblies folder. The script also produces a list of contigs extracted, which can be used later.

```bash
#Genome sized ACCs
$GETCIRCULARIZED $ASSEMBLY/assembly_info.txt $ASSEMBLY/assembly.fasta $ASSEMBLY/circularized_assemblies/circularized_to_analyze.txt $ASSEMBLY/circularized_assemblies 1000000

#Smaller ACCs
$GETCIRCULARIZED $ASSEMBLY/assembly_info.txt $ASSEMBLY/assembly.fasta $ASSEMBLY/circularized_assemblies/circularized_to_analyze.txt $ASSEMBLY/circularized_assemblies 0
```

---

## Polishing

The polishing pipeline (all the steps past this point) is available in Snakemake format. However, a step by step and some rationale is provided here as well.

Specify the contig number (Flye output) to polish.

```bash
CONTIG='1354'
```

### Minimap2

Align filtered reads:

```bash
minimap2 -ax map-ont $CIRC_CONTIGS/$CONTIG.fasta $FILTREADS > $MAPPED/$CONTIG.sam
```

### Gerenuq

Filter for high quality alignments to extract well mapped reads (by default, length > 1000, score > 1, length/score > 2, sequence identity > 0.9):

```bash
gerenuq -i $MAPPED/$CONTIG.sam -o $GQ_FILT/$CONTIG_gq.sam -t 8 -m 0.9

#convert gerenuq output (reads for polishing) to fasta (minipolish accepts fasta not sam)
samtools fasta $GQ_FILT/$CONTIG_gq.sam > $GQ_FILT/$CONTIG_gq.fasta
```

### Round 1: Minipolish

Racon + Medaka is a common method of polishing assemblies. We are using minipolish as it uses racon, but includes features to cleanly polish circular contigs.

Though minipolish expects miniasm assemblies, Flye is our chosen assembler as its explicitly designed for metagenomes. Flye assemblies can be made to work with minipolish with minor formatting modifications that don't directly involve modifying the sequence data.

1. A trailing "l" character is added to the end of the fasta header.
2. Convert the GFA format from GFA1 to GFA2
3. Skip initial is needed or minipolish won't run (this step - performing racon with reads used to build the contigs - requires non-standard data unique to miniasm)

```bash
#minipolish requires assembly graphs as input. Extract contig graphs from meta-assembly graph
$GFATOOLS view -l edge_361 -r 1 $ASSEMBLY/assembly_graph.gfa > $MINIPOL/$CONTIG.gfa

#convert to GFA2
gfak convert -S 2.0 $MINIPOL/$CONTIG.gfa > $MINIPOL/$CONTIG_gfa2.gfa

#miniasm formatting with a script (located in /scripts)
python3 $MINIASM_FORM -i $MINIPOL/$CONTIG_gfa2.gfa -o $MINIPOL/$CONTIG_asm.gfa
```

After the above modifications, minipolish can run.

```bash
#polish with gerenuq filtered reads
minipolish --skip_initial $GQ_FILT/$CONTIG_gq.fasta $MINIPOL/$CONTIG_asm.gfa > $MINIPOL/$CONTIG_polished.gfa
```

### Round 2: Medaka

```bash
#convert to fasta format
awk '/^S/{print ">"$2"\n"$3}' $MINIPOL/$CONTIG_polished.gfa | fold > $MINIPOL/$CONTIG_polished.fasta
#run medaka with gerenuq filtered reads
medaka_consensus -i $GQ_FILT/$CONTIG_gq.fasta -d $MINIPOL/$CONTIG_polished.fasta -o $MEDAKA/$CONTIG
```

Move all final assemblies to their own folder and name them appropriately

```bash
for fasta in $MEDAKA/*; do cp $fasta/consensus.fasta $FINAL_ASSEMBLY/${fasta##*/}.fasta; done
```

---

### QC: Coverage and Completeness

>We have our total reads mapped against our unpolished genomes already (from previous minimap2 steps), but now we want the coverage after polishing and filtering as well

### Polished and Filtered Coverage

Map gerenuq filtered reads against the polished assembly to get polished coverage

```bash
minimap2 -ax map-ont $FINAL_ASSEMBLY/$CONTIG.fasta $GQ_FILT/$CONTIG_gq.fasta > $MAP_POLISHED/$CONTIG.sam
```

Checking coverage over 1000 base windows with **Mosdepth**.

```bash
nohup samtools view -S -b $MAP_POLISHED/$CONTIG.sam -o $MAP_POLISHED/$CONTIG.bam &
samtools sort $MAP_POLISHED/$CONTIG.bam -o $MAP_POLISHED/$CONTIG_sort.bam
samtools index $MAP_POLISHED/$CONTIG_sort.bam

mosdepth --by 1000 $CONTIG $MAP_POLISHED/$CONTIG_sort.bam
```

Get final size (after polishing) and mean depth for all contigs in the folder from mosdepth output.

```bash
for file in $DEPTH_POLISH/*/*.mosdepth.summary.txt; do awk 'FNR==2 {printf $1 "\t"; printf $2 "\t" ;print $4}' ${file}; done > mosdepth_summary_stats.txt
```

### Completeness

**CheckM** completeness

```bash
nohup checkm lineage_wf -t 112 -x fasta $FINAL_ASSEMBLY $CHECKM &> $CHECKM/out.log &
```

### Running this on Secondary Assemblies

```bash
SECONDARY='Secondary_Pipeline'
│   └── {SAMPLES}
│       └── secondary_assembly
│           └── {BINS}
│                └── secondary_flye

SAMPLES=$SECONDARY/*
BINS=$SECONDARY/$SAMPLE/secondary_assembly/bin*
```

Handle secondary assembly data

```bash
#Iterates over and shows binned secondary assembly folders that have circularized contigs in them. Exclude contigs that have coverages less than 10 and are repetitive
for i in $ALL_BINS/assembly_info.txt; do awk '$3 < 10 && $4 == "Y" && $5 == "Y" {print $1}' $i; done | while read -r contig; do rm $contig; done

#This is just to handle the file structure of secondary assembly bins for running snakemake; renames assembly folder in each bin folder
GETCIRCULARIZED=/Volumes/data/gac_suncor/hsay/tools/extract-circularized-fasta.py
ASSEMBLY=3-assemblyhq

#rename folders to match the file structure snakemake workflow expects, go into each folder and extract circular contigs
for i in $BINS; do mv $i/secondary_flye $i/3-assemblyhq; mkdir $i/3-assemblyhq/circularized_assemblies; cd $i; $GETCIRCULARIZED $ASSEMBLY/assembly_info.txt $ASSEMBLY/assembly.fasta $ASSEMBLY/circularized_assemblies/circularized_to_analyze.txt $ASSEMBLY/circularized_assemblies 0; cd ..; done

#RUN SNAKEMAKE, set -d to bin folders. 
```
