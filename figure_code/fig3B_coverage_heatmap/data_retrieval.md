# Coverage Heatmap of 1% Bins

This contains code necessary for rerunning coverage calculations, then retrieves data output (coverage per base) from mosdepth results.

Rerun coverages if necessary:

```bash
#for every sample
for SAMP in *-*/INHERIT_phages; do \
 for a in $SAMP/*.fasta; do d=${a##*/}; \
  if [ ! -f $SAMP/coverage/gerenuq/${d%%.fasta}_gq.fasta ]; \
  then conda activate minimap2 && \
   minimap2 -ax map-ont $a /Volumes/backup/hsay/basecalls/${SAMP%%/*}/basecalled_filtered.fastq > $SAMP/coverage/minimap2/${d%%.fasta}.sam && conda activate gerenuq && \
   gerenuq -i $SAMP/coverage/minimap2/${d%%.fasta}.sam -o $SAMP/coverage/gerenuq/${d%%.fasta}_gq.sam -t 112 -m 0.9 -l 1000 && conda activate samtools && \
   samtools fasta $SAMP/coverage/gerenuq/${d%%.fasta}_gq.sam > $SAMP/coverage/gerenuq/${d%%.fasta}_gq.fasta && conda activate minimap2 && \
   minimap2 -ax map-ont $a $SAMP/coverage/gerenuq/${d%%.fasta}_gq.fasta > $SAMP/coverage/polished_coverage/${d%%.fasta}.sam && conda activate samtools && \
   samtools view -S -b $SAMP/coverage/polished_coverage/${d%%.fasta}.sam -o $SAMP/coverage/polished_coverage/${d%%.fasta}.bam && \
   samtools sort $SAMP/coverage/polished_coverage/${d%%.fasta}.bam -o $SAMP/coverage/polished_coverage/${d%%.fasta}_sort.bam && \
   samtools index $SAMP/coverage/polished_coverage/${d%%.fasta}_sort.bam && conda activate mosdepth && cd $SAMP/coverage/mosdepth && \
   mosdepth --by 1000 ${d%%.fasta} ../polished_coverage/${d%%.fasta}_sort.bam -t 112 && \
   cd /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids && \
   rm -r $SAMP/coverage/minimap2/${d%%.fasta}.sam $SAMP/coverage/gerenuq/${d%%.fasta}_gq.sam;\
  fi;\
  done;\
done

#unpolished coverage
for SAMP in *-*/INHERIT_phages; do \
 for a in $SAMP/*.fasta; do d=${a##*/}; \
  if [ ! -f $SAMP/coverage/mosdepth/${d%%.fasta}.per-base.bed ]; \
  then conda activate minimap2 && \
   minimap2 -ax map-ont $a /Volumes/backup/hsay/basecalls/${SAMP%%/*}/basecalled_filtered.fastq > $SAMP/coverage/unpolished_coverage/${d%%.fasta}.sam && conda activate gerenuq && \
   samtools view -S -b $SAMP/coverage/unpolished_coverage/${d%%.fasta}.sam -o $SAMP/coverage/unpolished_coverage/${d%%.fasta}.bam && \
   rm $SAMP/coverage/unpolished_coverage/${d%%.fasta}.sam
   samtools sort $SAMP/coverage/unpolished_coverage/${d%%.fasta}.bam -o $SAMP/coverage/unpolished_coverage/${d%%.fasta}_sort.bam && \
   rm $SAMP/coverage/unpolished_coverage/${d%%.fasta}.bam
   samtools index $SAMP/coverage/unpolished_coverage/${d%%.fasta}_sort.bam && conda activate mosdepth && cd $SAMP/coverage/mosdepth && \
   mosdepth --by 1000 ${d%%.fasta} ../unpolished_coverage/${d%%.fasta}_sort.bam -t 112 && \
   cd /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids && \
   rm -r $SAMP/coverage/unpolished_coverage/${d%%.fasta}.sam  $SAMP/coverage/unpolished_coverage/${d%%.fasta}_sort.bam;\
  fi;\
  done;\
done

#INHERIT phages folder
#we're keeping the polished_coverage folder
for d in *.fasta; do \
  if [ ! -f coverage/modepth/${d%%.fasta}.mosdepth.summary.txt ]; \
  then conda activate minimap2 && \
   minimap2 -ax map-ont $d /Volumes/backup/hsay/basecalls/10-T3_apr21-22/basecalled_filtered.fastq > coverage/minimap2/${d%%.fasta}.sam && conda activate gerenuq && \
   gerenuq -i coverage/minimap2/${d%%.fasta}.sam -o coverage/gerenuq/${d%%.fasta}_gq.sam -t 112 -m 0.9 -l 1000 && conda activate samtools && \
   samtools fasta coverage/gerenuq/${d%%.fasta}_gq.sam > coverage/gerenuq/${d%%.fasta}_gq.fasta && conda activate minimap2 && \
   minimap2 -ax map-ont $d coverage/gerenuq/${d%%.fasta}_gq.fasta > coverage/polished_coverage/${d%%.fasta}.sam && conda activate samtools && \
   samtools view -S -b coverage/polished_coverage/${d%%.fasta}.sam -o coverage/polished_coverage/${d%%.fasta}.bam && \
   samtools sort coverage/polished_coverage/${d%%.fasta}.bam -o coverage/polished_coverage/${d%%.fasta}_sort.bam && \
   samtools index coverage/polished_coverage/${d%%.fasta}_sort.bam && conda activate mosdepth && cd coverage/mosdepth && \
   mosdepth --by 1000 ${d%%.fasta} ../polished_coverage/${d%%.fasta}_sort.bam -t 112 && \
   cd ../../ && \
   rm -r coverage/minimap2/${d%%.fasta}.sam coverage/gerenuq/${d%%.fasta}_gq.sam;\
  fi;\
  done


```

Get and rename bed files into temp folder

```bash
for i in *-*; do \
 SNAME=${i##*/}
 for x in $i/INHERIT_phages/coverage/mosdepth/*.regions.bed*; do \
  cp $x ./temp/${SNAME}+${x##*/}; \
 done; \
done
```
