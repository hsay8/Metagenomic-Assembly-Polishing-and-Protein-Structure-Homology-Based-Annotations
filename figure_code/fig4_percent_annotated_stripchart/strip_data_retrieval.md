# **Compare number of actual gene products predicted vs hypotheticals**

The Foldseek and bakta output was used to generate the percent annotated stripcharts
Run the scripts within the INHERIT_phages folder of each sample

To compare the number of hypotheticals left between bakta and foldseek:

```bash
SAMPLE=

#count hypothetical products from foldseek and bakta output 
for i in *-*; do 
for x in $i/INHERIT_phages/*.fasta; do \
  f=${x##*/}; \
  #contig name
  echo -en "contig_${f%.fasta} \t"; \
  #number of cds predicted by bakta
  sed -ne  $'s/^.*CDSs: //p' $i/INHERIT_phages/bakta/${f%.fasta}.txt | tr -d '\n'; echo -ne "\t"; \
  #number of hypotheticals in bakta 
  sed -ne $'s/^.*hypotheticals: //p' $i/INHERIT_phages/bakta/${f%.fasta}.txt| tr -d '\n'; echo -ne "\t";\
  #number of hypotheticals identified via foldseek
  grep -c 'hypothetical' $i/INHERIT_phages/foldseek_hits/${f%.fasta}_besthits.tab;\
  done > $i/INHERIT_phages/${i##*/}_num_hypoths.txt ;\
done

```

Compare counts of GO/KEGG pathway annotations between bakta and foldseek:

```bash
for i in *-*; do \
for x in $i/INHERIT_phages/*.fasta; do \
  f=${x##*/}; \
  #contig name
  echo -en "contig_${f%.fasta} \t";\
  #number of cds
  sed -ne  $'s/^.*CDSs: //p' $i/INHERIT_phages/bakta/${f%.fasta}.txt | tr -d '\n'; echo -ne "\t"; \
  #number of KEGG and GO pathways annotated in bakta
  grep -c 'KEGG\|GO:' $i/INHERIT_phages/bakta/${f%.fasta}.gff3 | tr -d '\n'; echo -ne "\t"; \
  #number of GO pathways annotated in foldseek
  grep -c 'GO:' $i/INHERIT_phages/foldseek_hits_to_anno/${f%.fasta}_besthits_annot.txt| tr -d '\n'; echo -ne "\t"; \
  #count number of KEGG pathways annotated in foldseek
  awk -v num="$n" -F"\t" '$10{c++}END{print c-1+$n}' $i/INHERIT_phages/foldseek_hits_to_anno/${f%.fasta}_besthits_annot.txt; \
  done > $i/INHERIT_phages/${i##*/}_num_path_annos.txt ;\
done
```

Rerun coverage calculations for phages

```bash
nohup bash -c "for x in /Volumes/backup/hsay/basecalls/*; do \
  for y in ./${x##*basecalls/}/INHERIT_phages/*.fasta; do \
  i=${y##*/}; \
    minimap2 -ax map-ont $y $x/basecalled_filtered.fastq > ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}.sam; \
      samtools view -S -b ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}.sam -o ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}.bam; \
      samtools sort ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}.bam -o ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}_sort.bam; \
      samtools index ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}_sort.bam; \
    mosdepth --by 1000 ./${x##*basecalls/}/INHERIT_phages/coverage/mosdepth/${i%.fasta} ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}_sort.bam; \
      rm ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}.sam ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}.bam ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}_sort.bam ./${x##*basecalls/}/INHERIT_phages/coverage/minimap2/${i%.fasta}_sort.bam.bai;\
   done;\
   done" & 

```