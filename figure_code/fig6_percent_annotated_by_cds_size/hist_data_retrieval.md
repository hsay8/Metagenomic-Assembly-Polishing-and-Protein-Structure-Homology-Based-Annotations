# Annotation likelihood by size - compare likelihood of annotation vs CDS size by bakta and foldseek

```bash
#Get locus tags and their sizes
for f in *-*; do \
  for i in $f/INHERIT_phages/bakta/*.tsv; do \
  #skip over hypotheticals.tsv's \
    if [[ "$i" != *"hypotheticals"* ]]; then
    #print locus tag and CDS size
    awk -F "\t" 'NR>3 {printf $6; printf "\t"; print $4-$3}' $i ; \
    fi ; 
  done > $f/INHERIT_phages/${f##*/}_CDS_by_size.txt; \
done

#Get locus tags and their foldseek annotations
for f in *-*; do \
  for i in $f/INHERIT_phages/foldseek_hits_to_anno/*.txt; do \
  #print locus tag, GO, KEGG gene
    awk -F"\t" '$1!=" " && NR>1{print substr($1,1,12),"\t", "Y"}' $i; \
    done > $f/INHERIT_phages//${f##/}_fsannot_by_tag.txt ;\
done 

#Get locus tags and their bakta annotations - except we're getting hypotheticals
#I dont wan't to deal with negative matching with awk. It'll be easier to get loci with no annotations to find ones with annotations in R.
for f in *-*; do \
  for i in $f/INHERIT_phages/bakta/*.tsv; do \
  #skip over hypotheticals.tsv's \
    if [[ "$i" != *"hypotheticals"* ]]; then
    #print locus tag and CDS size
    awk -F "\t" 'NR>3 && $8=="hypothetical protein"{print $6, "\t", "N"}' $i ; \
    fi ; 
  done > $f/INHERIT_phages/${f##*/}_bkannot_by_tag.txt; \
done
```
#transfer to local directory (windows)
```cmd
scp -r hsay@129.100.236.35:\Volumes\data\gac_suncor\hsay\gac_samples\allvall_plasmids\*-*\INHERIT_phages\*by_tag* C:\Users\hh_sa\OneDrive\GitHub\Suncor-Project-\figure_code\percent_annotated_by_cds_size\data
```

Get coverage by CDS size
I'll get the predicted CDS in a contig, and the per base coverage (mosdepth) files separately to merge in R. 
```bash
for f in *-*; do \
  for i in $f/INHERIT_phages/bakta/*.tsv; do \
  #skip over hypotheticals.tsv's \
    if [[ "$i" != *"hypotheticals"* ]]; then
    #print locus tag and CDS start and stop 
    awk -F "\t" 'NR>3 && $8=="hypothetical protein"{print $6, "\t", $3, "\t", $4}' $i ; \
    fi ; 
  done > $f/INHERIT_phages/${f##*/}_bkannot_ranges.txt; \
done

```