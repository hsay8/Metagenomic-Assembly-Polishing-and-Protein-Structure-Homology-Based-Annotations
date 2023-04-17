## **KEGG Pathway Completeness**
We'll compare the pathway completeness of bakta and foldseek by first getting the KO's from the output. Then we'll attribute the KO's to a pathway. Then we'll retrieve the complete list of KO's for pathways found. The KO's we have from the output that exist in that list will tell us the completeness of the pathway. 

Get KEGG KO's from **bakta** output
```bash
for f in *-*; do \
  for i in $f/INHERIT_phages/bakta/*.tsv; do \
  #will look for KEGG: or GO: matches, will print whatever's between KEGG: or GO: and a comma
    grep -oP '(?<=KEGG:).*?(?=,)' $i; \
    done > $f/INHERIT_phages/${f##/}_all_phages_bakta_KO.txt ;\
  done 
```

Retrieve KEGG pathways data from database using KO (by bakta)
```bash
for i in *-*; do \
  for f in $i/INHERIT_phages/${i##/}_all_phages_bakta_KO.txt; do \
    while read -r x; do curl -s "https://rest.kegg.jp/link/pathway/${x## }"; done < $f; \
    done > $i/INHERIT_phages/${i##/}_all_phages_bakta_KOs_to_path.txt ;\
done
```

Get full list of KO's from pathways found
```bash
for i in *-*; do \
  for f in $i/INHERIT_phages/${i##/}_all_phages_bakta_KOs_to_path.txt; do \
    while IFS=$'\t' read -r x y; do curl -s "https://rest.kegg.jp/get/${y## }/"; done < $f; \
    done > $i/INHERIT_phages/${i##/}_all_phages_bakta_pathways_all_KOs.txt ; \
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
Retrieve data from KEGG database using KEGG gene identifiers provided by foldseek 
```bash
#Get KO's from KEGG gene identifiers and pathways from KO
for i in *-*; do \
  for f in $i/INHERIT_phages/${i##/}_all_phages_foldseek_function.txt; do \
    awk -F "\t" '$3!=" " {gsub(";","+"); print $3}' $f | while IFS=$'\t' read -r x; do curl -s "https://rest.kegg.jp/link/ko/${x## }"; done | sort -u; \
    done > $i/INHERIT_phages/${i##/}_all_phages_foldseek_KOs.txt; \
done

#Get pathways from KO's
for i in *-*; do \
  for f in $i/INHERIT_phages/${i##/}_all_phages_kegg-gene_KOs.txt; do \
    while IFS=$'\t' read -r x y; do curl -s "https://rest.kegg.jp/link/pathway/${y## }" ; done < $f; \
    done > $i/INHERIT_phages/${i##/}_all_phages_foldseek_KOs_to_path.txt ;\
done

#Get full list of KO's from pathways
for i in *-*; do \
  for f in $i/INHERIT_phages/${i##/}_all_phages_foldseek_KOs_to_path.txt; do \
    while IFS=$'\t' read -r x y; do curl -s "https://rest.kegg.jp/link/ko/${y## }"; done < $f | sort -u; \
    done > $i/INHERIT_phages/${i##/}_all_phages_foldseek_pathways_all_KOs.txt ;\
done

```

KEGG pathway modules contain KO's that represent a function in the pathway