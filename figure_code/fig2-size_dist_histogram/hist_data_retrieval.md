
# **HISTOGRAM**

Get the lengths of each contig into a list using sequence stats

```bash
for f in *-*; do cat $f/initial_polished/*.fasta $f/secondary_polished/*.fasta > ${f}_concat.fasta; done

for f in *_concat.fasta; do /Volumes/data/gac_suncor/hsay/tools/sequence-stats-1.0 -c $f > ${f%%_concat*}_size_dist.txt; done
```

For phages (size after polishing, uses mosdepth output)

```bash
for f in *-*; do for i in $f/INHERIT_phages/coverage/mosdepth/*.summary.txt; do awk 'FNR == 2 {printf $1; printf "\t"; print $2}' $i ; done > $f/INHERIT_phages/${f##*/*-}_size_dist.tab;  done
```

USING BIOAWK (OLD)

Get the lengths of each contig into a list

```bash
##actually, the naming doesn't even really matter
#concat circles
for f in *-*; do x=${f#*-}; for c in $f/initial_polished/*.fasta; do cat $c >> $x.fasta; done; for c in $f/secondary_polished/*.fasta; do cat $c >> $x.fasta; done; done

#label fasta headers with the sample its from using filenames as the sample name
for x in *.fasta; do sed "s/>.*/&;${x%.fasta}/" $x > ${x%%.fasta}_2.fasta; done

#using bioawk for all circles
for c in *.fasta; do bioawk -c fastx '{ print $name, length($seq) }' < $c >> ./all_contigs_plasmid_sizes.txt; done

#for all phages
for c in *-*/INHERIT_phages/*.fasta; do bioawk -c fastx '{ print $name, length($seq) }' < $c >> ./all_phages_sizes.txt; done
```

Get recurring plasmids from blast results

```bash
#for all recurring plasmid (no phages)

for i in /Volumes/data/gac_suncor/hsay/gac_samples/allvall_plasmids/BLAST/recurring_plasmids/*/; do for c in $i/final_assemblies/*.fasta; do bioawk -c fastx '{ print $name, length($seq) }' < $c >> ./recurring_plasmid_sizes.txt; done ;done

##for recurring phages
#this script just gets all the phages from the recurring plasmids using INHERIT's output
for i in *-*/; do while read phage type score; \do x=${i#*-}; mv ./BLAST/recurring_plasmids/$x/final_assemblies/$phage.fasta ./BLAST/recurring_plasmids/$x/final_assemblies/phages/$phage.fasta; done < $i/INHERIT_predicted_phages.txt; done
#gets recurring phages
for i in /Volumes/data/gac_suncor/hsay/gac_samples/allvall_plasmids/BLAST/recurring_plasmids/*/; do for c in $i/final_assemblies/phages/*.fasta; do bioawk -c fastx '{ print $name, length($seq) }' < $c >> ./recurring_phages_sizes.txt; done ;done
```
