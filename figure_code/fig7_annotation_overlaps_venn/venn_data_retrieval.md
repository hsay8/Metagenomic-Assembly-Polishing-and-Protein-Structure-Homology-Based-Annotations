## Result Matching Between Foldseek and Bakta
Here we compare the bakta and foldseek annotations based on locus tags. This will be used to compare the two outputs to see if they agree on the annotation of a predicted CDS. 

Tables will be generated for bakta (containing the tag, and accession of hits), then for foldseek.
The merge function in R will be helpful in this case to get proteins which they BOTH have predictions for, and we can compare to see if they agree. 


### Comparing PFAMS
Get PFAMS from **bakta** output 
```bash
for f in *-*; do \
  for i in $f/INHERIT_phages/bakta/*.tsv; do \
  #skip over hypotheticals.tsv's
    if [[ "$i" != *"hypotheticals"* ]]; then
      awk -F "\t" -v OFS=, '$9 ~ /PFAM:/ {print $6, "/t", $9}' $i | \
      #look for PFAM, strip any whitespaces, print locus tag and PFAM.
      awk -F, '{for(i=1;i<=NF;i++){if($i~/PFAM:/){\
        gsub("PFAM:", "", $i);\
        print $1, "\t", $i }\
        }}' ;\
      fi;
    done > $f/INHERIT_phages/${f##/}_all_phages_bakta_pfam.txt ; \
  done 
```

Get PFAMS from **foldseek** outputs
```bash
for f in *-*; do \
  for i in $f/INHERIT_phages/foldseek_hits_to_anno/*.txt; do \
  #print locus tag, PFAM | remove rows with empty pfam | 
    awk -F"\t" '$1!=" " && NR>1{\
      print substr($1,1,12),"\t", substr($11, 1, length($11)-1)}' $i | awk '$2!=""' | sort -u ; \
    done > $f/INHERIT_phages/${f##/}_all_phages_foldseek_pfam.txt ;\
done 
```


### Comparing protein product predictions
Retrieves bakta output containing a UniRef(50 or 90) match, keeping the 90 hit over the 50 if it exists.
```bash
for f in *-*/INHERIT_phages/bakta/; do \
for i in $f/*.tsv; do \
#skip over hypotheticals.tsv's
if [[ "$i" != *"hypotheticals"* ]]; then
#First print locus tag and annotations that have a UniRef hit. Since annotations that we don't need for this are jammed into one column (and separated by "," i.e "KEGG:K07097, SO:0001217, UniRef:UniRef50_A0A7X4AT10, UniRef:UniRef90_A0A6N0ZF79)' change the delimiter to a comma and separate later. 
awk -F "\t" -v OFS=, '$9 ~ /UniRef:/ {print $6, $9}' $i | \
#Have awk recognize the comma as a delimiter to break up the long list of annotations (from the example above) from a single column/cell into individual columns, then print only the whole match (UniRef:) and corresponding locus tag (col 1). Gsub to remove prefixes.
awk -F, '{for(i=1;i<=NF;i++){if($i~/UniRef:/){\
  gsub("UniRef:UniRef50_", "", $i);\
   gsub("UniRef:UniRef90_", "", $i);\
    gsub("UniRef:UniRef100_", "", $i);\
     print $1, $i }\
   }}' | \
#sort by accession column in reverse, then get uniques (by column 1). In this order, uniref 90 hits will be kept instead of 50 if it exists 
sort -r -k2,2 | sort -u -k 1,1; \
fi; done > ${f%/bakta/}/${SAMPLE}_bakta_hits_accessions.txt;\
done
```

From the Foldseek output (Alphafold/Uniprot), get the accession of hits.
```bash
for f in *-*/INHERIT_phages/foldseek_hits_to_anno/; do \
for i in $f/*.txt; do \
#print locus tag in col 1, then accession in col 2
awk '{print substr($1,1,12), $2}' $i; \
done > ${f%/foldseek_hits_to_anno/}/${SAMPLE}_foldseek_hits_accessions.txt;\
done
```

### Comparing GO terms and KEGG KO's
This data is used for generating the venn diagram, which sees if bakta annotations and annotations obtained via foldseek is in agreement.

Get KO's, GO terms from **bakta** output
```bash
for f in *-*; do \
  for i in $f/INHERIT_phages/bakta/*.tsv; do \
  #will look for KEGG: or GO: matches, will print whatever's between KEGG: or GO: and a comma
    grep -oP '(?<=KEGG:).*?(?=,)' $i; \
    grep -oP '(?<=GO:).*?(?=,)' $i;\
    done > $f/INHERIT_phages/${f##/}_all_phages_bakta_function.txt ;\
  done 
```

Get KO's, GO terms from **foldseek hits** (from which GO and KEGG GENES have been retrieved from uniprot)
KEGG gene names will be converted to KO's in the next step
```bash
for f in *-*; do \
  for i in $f/INHERIT_phages/foldseek_hits_to_anno/*.txt; do \
  #print locus tag, GO, KEGG gene
    awk -F"\t" '$1!=" "&&NR>1{print substr($1,1,12),"\t", $9, "\t", substr($10, 1, length($10)-1)}' $i | awk '$2!=""&&$3!=""' | sort -u ; \
    done > $f/INHERIT_phages//${f##/}_all_phages_foldseek_function.txt ;\
done 
```
Convert kegg gene identifiers to KO's
```bash
for i in *-*; do \
  for f in $i/INHERIT_phages/${i##/}_all_phages_foldseek_function.txt; do \
    awk -F "\t" '$3!=" " {gsub(";","+"); print $3}' $f | while IFS=$'\t' read -r x; do curl -s "https://rest.kegg.jp/link/ko/${x## }"; done | sort -u; \
    done > $i/INHERIT_phages/${i##/}_all_phages_kegg-gene_KOs.txt ;\
done
```

