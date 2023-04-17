## GO Terms Pie Chart

Retrieves GO terms from foldseek annotations.

Get go terms into one file

```bash
for i in *-*/INHERIT_phages/foldseek_hits_to_anno/*_besthits_annot.txt; do awk -F "\t" 'NR>1 && NF {print $8}' $i | grep -v '^$'; done > goterms_total.txt
```

R doesn't like single quotes. It messes up read.table, so remove them.

```bash
sed -i "s/'//g" goterms_total.txt
```
