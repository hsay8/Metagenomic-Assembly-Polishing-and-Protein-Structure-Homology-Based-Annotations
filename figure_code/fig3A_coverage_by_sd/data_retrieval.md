# Rerun coverages if necessary

```bash

#phages
SAMPLE="2-T4_may14-21"

for i in $SAMPLE/*_polished; do for c in *.fasta; \
 do d=${c##*/}; \
  if [ ! -f modepth/${d%%.fasta}.mosdepth.summary.txt ]; \
  then conda activate minimap2 && \
   minimap2 -ax map-ont $c /Volumes/backup/hsay/basecalls/${i%%/*}/basecalled_filtered.fastq > minimap/${d%%.fasta}.sam && conda activate gerenuq && \
   gerenuq -i minimap/${d%%.fasta}.sam -o gerenuq/${d%%.fasta}_gq.sam -t 112 -m 0.9 && conda activate samtools && \
   samtools fasta gerenuq/${d%%.fasta}_gq.sam > gerenuq/${d%%.fasta}_gq.fasta && conda activate minimap2 && \
   minimap2 -ax map-ont $c gerenuq/${d%%.fasta}_gq.fasta > polished_coverage/${d%%.fasta}.sam && conda activate samtools && \
   samtools view -S -b polished_coverage/${d%%.fasta}.sam -o polished_coverage/${d%%.fasta}.bam && \
   samtools sort polished_coverage/${d%%.fasta}.bam -o polished_coverage/${d%%.fasta}_sort.bam && \
   samtools index polished_coverage/${d%%.fasta}_sort.bam && conda activate mosdepth && cd mosdepth && \
   mosdepth --by 1000 ${d%%.fasta} /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids/polished_coverage/${d%%.fasta}_sort.bam -t 112 && \
   rm -r /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids/minimap/${d%%.fasta}.sam /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids/gerenuq/${d%%.fasta}_gq.sam /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids/polished_coverage/${d%%.fasta}.sam && \
   cd /Volumes/data/gac_/hsay/gac_samples/allvall_plasmids/;\
  fi;\
  done;\ 
done 
```

Get the completeness of recovered genomes - these will basically be the benchmark for good quality circles when looking at phages

```bash
#find assemblies that are likely genomes: aka above 500 kilobytes fasta. move to a temp folder
find *-*/final_assemblies -name "*.fasta" -size +500k -exec echo {} \; | while read i; do cp $i ./temp/${i%%/*}+${i##*/}; done


#run checkm on temp folder 
nohup checkm lineage_wf -t 112 -x fasta ./temp ./temp/checkm &> ./temp/checkm/completeness.log &
cd checkm
checkm qa ./lineage.ms ./ -t 112 -f completeness.tab

#subset for completeness over 90, contam less than 10 
head -n -1 ./checkm/completeness.tab | awk -F " " 'NR>3 && $13 > 90 && $14 < 10 {print $1,"\t", $13, "\t", $14}' | while read x y z; do \
#sample
i=${x%%+*};\
#contig
d=$x.fasta
#run coverage
conda activate minimap2 && \
 minimap2 -ax map-ont $d /Volumes/backup/hsay/basecalls/${i%%/*}/basecalled_filtered.fastq > minimap/${d%%.fasta}.sam && conda activate gerenuq && \
 gerenuq -i minimap/${d%%.fasta}.sam -o gerenuq/${d%%.fasta}_gq.sam -t 112 -m 0.9 && conda activate samtools && \
 samtools fasta gerenuq/${d%%.fasta}_gq.sam > gerenuq/${d%%.fasta}_gq.fasta && conda activate minimap2 && \
 minimap2 -ax map-ont $d gerenuq/${d%%.fasta}_gq.fasta > polished_coverage/${d%%.fasta}.sam && conda activate samtools && \
 samtools view -S -b polished_coverage/${d%%.fasta}.sam -o polished_coverage/${d%%.fasta}.bam && \
 samtools sort polished_coverage/${d%%.fasta}.bam -o polished_coverage/${d%%.fasta}_sort.bam && \
 samtools index polished_coverage/${d%%.fasta}_sort.bam && conda activate mosdepth && cd mosdepth && \
 mosdepth --by 1000 ${d%%.fasta} /Volumes/data/gac_/hsay/temp/polished_coverage/${d%%.fasta}_sort.bam -t 112 && \
 rm -r /Volumes/data/gac_/hsay/temp/minimap/${d%%.fasta}.sam /Volumes/data/gac_/hsay/temp/gerenuq/${d%%.fasta}_gq.sam /Volumes/data/gac_/hsay/temp/polished_coverage/${d%%.fasta}.sam && \
 cd /Volumes/data/gac_/hsay/temp;\
done

```
