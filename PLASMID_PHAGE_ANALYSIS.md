# Plasmid and Phage Analysis

- [Plasmid and Phage Analysis](#plasmid-and-phage-analysis)
  - [All vs. All BLAST of ACCs \< 1mb](#all-vs-all-blast-of-accs--1mb)
  - [Annotation - By Sequence Homology](#annotation---by-sequence-homology)
    - [Bakta](#bakta)
  - [Phage Identification and Classification](#phage-identification-and-classification)
  - [Annotation - By Protein Structure](#annotation---by-protein-structure)
    - [Colabfold - Hypothetical Protein Structure and Complex Prediction](#colabfold---hypothetical-protein-structure-and-complex-prediction)
    - [Structure Homology Search](#structure-homology-search)
  - [Foldseek Phage Structural Protein Clustering](#foldseek-phage-structural-protein-clustering)

**Directory Tree:**

The naming of each assembly (assembly #) is preserved from what Flye had assigned them (edge_#). Secondary assebmlies additionally contain a prefix including the bin it comes from.

```bash
WKDIR='Small_ACCs'
├── {SAMPLES} #Each folder contains assembled and polished ACCs from a given sample
│   ├── BLAST
│   ├── INHERIT_phages
│   │   ├── bakta
│   │   ├── colabfold_predictions
│   │   ├── foldseek_hits
│   │   ├── foldseek_hits_to_anno
│   │   └── {CONTIGS}.fasta
│   ├── initial_polished #assembled with just flye
│   │   ├── bakta
│   │   ├── has_phage_protein
│   │   ├── contains_ribosomal
│   │   ├── has_phages_INHERIT_results.txt #contains list of scores for ACC's with predicted phage proteins (bakta)
│   │   └── {CONTIGS}.fasta
│   └── secondary_polished #assembled with secondary pipeline
│       ├── bakta
│       ├── has_phage_protein
│       ├── contains_ribosomal
│       ├── contains_ribosomal
│       ├── has_phages_INHERIT_results.txt
│       └── {CONTIGS}.fasta

#Sample directories
SAMPLES=$WKDIR/*

#Polished contigs
ALL_INIT_CONTS=$WKDIR/$SAMPLE/initial_polished/*.fasta
ALL_SEC_CONTS=$WKDIR/$SAMPLE/secondary_polished/*.fasta

#PHAGES
BK_PRED_PHAGES=$WKDIR/$SAMPLE/
ALL_PHAGES=$WKDIR/$SAMPLE/*.fasta
PHAGE_FS=$WKDIR/$SAMPLE/INHERIT_phages/

```

## All vs. All BLAST of ACCs < 1mb

For finding singleton/recurring assemblies, create BLAST database of ACCs across GAC samples. ACCs that have high sequence identity between samples are considered the same (recurring).

```bash
#Before creating the database, add a suffix to headers of secondary assemblies (named by {bin}_{contig number}.fasta) to indicate what bins they come from. Currently with the way that the secondary assembly works, each secondary assembly will likely have the same edge number, so this is needed to distinguish them after using doing the BLAST.
for x in $ALL_SEC_CONTS; do y=${x##*basecalls/}/INHERIT_phages; sed "s/>.*/&_${y%%_*}/" $x > ${x%%.fasta}_headers.fasta; done

#concat fastas (initial and secondary) from each sample into sample specific multi fastas
for s in $SAMPLES; do cat $s/initial_polished/*.fasta > ./${s}_concat.fasta; cat $s/secondary_polished/*_headers.fasta >> ./${s}_concat.fasta; mv *.fasta $WKDIR/BLAST/; done

#in the BLAST folder, add a suffix (of the fasta filename, which will be the sample name) to headers, use ";" as a separator.
for FASTAS in *.fasta; do y=${x#*-}; sed "s/>.*/&;${y%_concat.fasta}/" $x > ${x%%.fasta}_2.fasta; done

#then make blast db
cat *_2.fasta > database.fasta
makeblastdb -in database.fasta -out gac_plasmids_BLASTdatabase -dbtype nucl
```

ALLvALL BLAST

```bash
##PLASMIDS
blastn -db gac_plasmids_BLASTdatabase -query ./database.fasta -outfmt "7 qacc qlen sacc slen evalue qstart qend sstart send qcovs nident" -num_threads 112 -out results.tab

#filter blast results to see recurring plasmids: percent identity >= 98%, query coverage within 10% of query length, query length is at least 90% of subject length, remove same plasmid matching and repeats
sed '/#/d' ./results.tab | awk '{ if (($10 >= 98) && (($11 / $2) > 0.90) && (($2 / $4 < 1.10) && ($2 / $4 >0.90) && ($1 != $3)) ) { print } }' > recurring_plasmids.tab

#Less stringent (90% identity, query not the same as subject)
sed '/#/d' ./results.tab | awk '{ if (($10 > 90) && (($11 / $2) > 0.80) && ($1 != $3)) { print } }' > results_90iden.tab
```

## Annotation - By Sequence Homology

Predict CDS with Bakta, then pass them to Colabfold to predict the structure. Then use Foldseek to query the best predicted structures for each CDS.

### Bakta

Using bakta

```bash
#From working directory, run bakta on all samples in intial/secondary polished folders.
nohup bash -c 'for SAMPLE in *; do for c in $SAMPLE/initial_polished/*.fasta; do bakta --db /Volumes/data/gac_suncor/hsay/gac_samples/bakta_db/db $c -o $SAMPLE/initial_polished/bakta -t 112 --complete; done; done' &

nohup bash -c 'for SAMPLE in *; do for c in $SAMPLE/secondary_polished/*.fasta; do bakta --db /Volumes/data/gac_suncor/hsay/gac_samples/bakta_db/db $c -o $SAMPLE/secondary_polished/bakta -t 112 --complete; done; done' &

#Get list of files containing ribosomal annotations from bakta output
for ANNO in */*_polished/bakta; do \
  grep --include=\*.gff3 -rwHl $ANNO -e 'ribosomal' > $ANNO/contigs_to_exclude.txt ;\
done 

#From working directory, move fastas of contigs listed in contigs_to_exclude.txt
for i in */*_polished/; do \
while read c; do
  c="${c##*/}";
  mv $i/${c%%.gff3}.fasta $i/contains_ribosomal/;
  done < $i/bakta/contigs_to_exclude.txt ;\
done 

```

## Phage Identification and Classification

With the bakta results, get list of contigs that contain phage proteins and move phages into the has_phage_protein folder in the initial and secondary assemblies folder.

```bash
#From working directory 
for s in */*_polished; do for i in $s/bakta/*.gbff; do awk '/phage/ {print FILENAME}' $i; done | while read -r contig; do x=${contig##*/}; mv $s/${x%%.gbff}.fasta $s/has_phage_protein/; done; done
```

Concat phage fastas into a multifasta for each sample bins to run INHERIT

```bash
#concat phages
for s in *-*; do cat $s/initial_polished/has_phage_protein/*.fasta $s/secondary_polished/has_phage_protein/*_headers.fasta > $s/all_w_phage_prot.fasta; done
```

INHERIT (phage identification: deep learning + database method)
Use pre-trained:

```bash
python3 IHT_predict.py --sequence $WKDIR/*_polished/all_w_phage_prot.fasta --withpretrain True --model INHERIT.pt --out $WKDIR/*/has_phages_INHERIT_results.txt
```

Copy all predicted phages from INHERIT OUTPUT into a new folder

```bash
#from working directory
for i in *-*; do grep 'Phage' $i/has_phages_INHERIT_results.txt > $i/INHERIT_predicted_phages.txt; done

for s in *-*; do while read contig type score; do c=${contig#edge_}; e=${c##*_}; p=${c%%l*}; cp $s/initial_polished/has_phage_protein/${c%l*}.fasta $s/INHERIT_phages/; cp $s/secondary_polished/has_phage_protein/${contig}.fasta $s/INHERIT_phages/; done < $s/INHERIT_predicted_phages.txt; done

for i in INHERIT_phages/*bin*.fasta; do x=${i##*/}; mv ../secondary_polished/bakta/${x%.fasta}.faa INHERIT_phages/bakta/; done 
```

## Annotation - By Protein Structure

### Colabfold - Hypothetical Protein Structure and Complex Prediction

```bash
#run for each sample folder in WKDIR
for f in */; do for i in $f/INHERIT_phages/bakta/*.faa; do x=${i##*/}; \
  /programs/colabfold_batch/bin/colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax $i $f/INHERIT_phages/colabfold_predictions/${x%.faa}; done &> $f/INHERIT_phages/colabfold_predictions/colab.log & done
```

### Structure Homology Search

FOLDSEEK (run per sample basis - within each sample folder)

```bash
mkdir foldseek_hits foldseek_hits_to_anno

#Runs foldseek and puts the best hits into a new tab file. Since the results are already sorted by best hit (sorted by e-val), we can use "sort -u" which keeps only first occurance of an entry of a column (which would be the best hit)
#Run from within sample folder
nohup bash -c 'for f in colabfold_predictions/*; do \
  foldseek easy-search $f/*unrelaxed_rank_1*.pdb /Volumes/data2/database/afdb ./foldseek_hits/${f##*/}.tab ./foldseek_hits/tmp --format-mode 0 --format-output "query,target,evalue,pident,fident,bits,qcov,tcov,alntmscore" -e 0.001; sort -u -k1,1 ./foldseek_hits/${f##*/}.tab > ./foldseek_hits/${f##*/}_besthits.tab ; \
  done' &> ./foldseek_hits/foldseek.log &
```

Getting functional annotations from foldseek hits (runs for all samples, outermost directory)

```bash
#return fields from UNIPROT
for i in *-*/INHERIT_phages/foldseek_hits/; do
  for f in $i/*_besthits.tab; do \
    #print heaters and initialize file
    file=${f##*/}; 
    echo -e "Query\tTarget Accession\tProtein Name\tFunction [CC]\tPathway [CC]\tAnnotation Score\tProtein families\tGene Ontology (GO)\tGene Ontology IDs\tKEGG\tPFAM" > ${i%foldseek_hits/}/foldseek_hits_to_anno/${file%.tab}_annot.txt; \
    #read foldseek best hits and query them in uniprot
    while read query target evalue pident fident bits qcov tcov alntmscore; do acc=${target#AF-}; \
    #print query name first 
        echo -en "${query%%_unrelaxed*} \t" >> ${i%foldseek_hits/}/foldseek_hits_to_anno/${file%.tab}_annot.txt; \
        #return fields 
        curl "https://rest.uniprot.org/uniprotkb/search?query=accession:${acc%%-*}&format=tsv&fields=accession,protein_name,cc_function,cc_pathway,annotation_score,protein_families,go,go_id,xref_kegg,xref_pfam" | sed -n 2p >> ${i%foldseek_hits/}/foldseek_hits_to_anno/${file%.tab}_annot.txt ; done < $f; \
  done; 
done
```

## Foldseek Phage Structural Protein Clustering

To try and relate phage, we take any Foldseek annotated phage head/tail/capsid proteins and cluster them with predicted structures across all samples.

After blasting, I copied the recurring phage into its own work directory.

```bash
WKDIR=./recurring_phages
├── recurring_phages
│   └── data
│      ├── {recurring_phages}.fasta
│      ├── structures
│      │ └── {contig_folders}
│      │  └── {proteins_colabfold_output_folders}
│      │   └── {protein}.pdb
│      ├── annotations
│       └── foldseek
│      ├── {db files}
│      └── aln
│       └── {alignment_files}
```

Get protein names representatives/members containing the keyword phage

```bash
#clusters all the proteins predicted
pdbs=./data/structures/*/*.pdb
db=./data/structures/db
clst=./data/clustering/clst
aln=./data/structures/aln/aln

##replace names in tsv results with foldseek's best hit + hits accessions if they exist
#make a copy to edit
cp ./data/clustering/clu.tsv ./data/clustering/clu_fsnames.tsv
#first get bakta name, hit name, accession from foldseek annotations. Then use sed to replace bakta names in our cluster results were foldseek names exist
for i in ./data/annotations/foldseek/*_annot.txt; do \
  awk -F "\t" 'NR>1{print $1, $2, $3}' $i | while read bk_prot acc hit_name; do \
    hit_name=$(echo $hit_name | tr -s '[:blank:]' '_') ;\
    sed -i --expression "s@${bk_prot}[^[:space:]]*.pdb@${acc}+${hit_name}@g" ./data/clustering/clu_fsnames.tsv ;\
done ;\
done

#pull entries containing the keyword phage
grep -i -E 'phage|tail|head|capsid|baseplate' ./data/clustering/clu_fsnames.tsv > ./data/clustering/clu_phage_prot.tsv

#Reverse match to get protein sequences from bakta. 
#Take the cluster representative names (which is the accession), match it to foldseek results accession column, pull bakta name, use bakta names to get files
cat ./data/clustering/clu_phage_prot.tsv | \
while IFS="+" read rep mem; do \
  #some representatives don't have FS annotations and wont have correct names, in that case just take the bakta name directly by grabbing the first 12 letters (bakta naming is only 12 letters, accessions are always less than 12)
  rep=${rep:0:12} ;\
  #find foldseek results files that contain the representative accession, cp those files into a new folder.
  find ./data/annotations/foldseek/ -type f -iname "*_annot.txt" -exec grep -ln "$rep" {} + | while read x; do x=${x%%_besthits_annot.txt}; cp /Volumes/data/gac_suncor/hsay/gac_samples/allvall_plasmids/*-*/INHERIT_phages/bakta/${x##*/}.faa /Volumes/data/gac_suncor/hsay/gac_samples/allvall_plasmids/recurring_phages/data/annotations/bakta; done ;\
done

#get bakta names back of phage proteins reps only so we can pull the sequences later
cat ./data/clustering/clu_phage_prot.tsv | awk '{print $1}' | sort -u | \
while IFS="+" read acc ignore; do \
  acc=${acc:0:12} ;\
  find ./data/annotations/foldseek/ -type f -iname "*_annot.txt" -exec grep -ln "$acc" {} + | while read x; do awk -F "\t" -v pat="$acc" '$2 ~ pat && NR>1 {print $1}' $x;\
done; done | sort -u > ./data/clustering/reps_phage_bkprot.tsv

##also get bakta names of the members to phage reps
cat ./data/clustering/clu_phage_prot.tsv | awk '{print $2}' |\
while IFS="+" read acc ignore; do \
  acc=${acc:0:12} ;\
  find ./data/annotations/foldseek/ -type f -iname "*_annot.txt" -exec grep -ln "$rep" {} + | while read x; do awk -F "\t" -v pat="$acc" '$2 ~ pat && NR>1 {print $1}' $x;\
done; done | sort -u > ./data/clustering/mems_phage_bkprot.tsv
```

Get AA sequences, going to rename the bakta named headers to foldseek hit accession

```bash
#Make sure we have our bakta annotations in our new working directory
for i in *.fasta; do find ../../* -type f -iname "${i%%.fasta}.faa" -exec cp -t ./annotations/bakta/ {} +; done

#Get AA sequences of all reps and members after clustering
cat ./data/clustering/clu.tsv |\
while read rep mem; do \
  rep=${rep:0:12} ;\
  mem=${mem:0:12} ;\
  find ./data/annotations/bakta/ -type f -iname "*.faa" -exec grep -ln "$rep" {} + | \
  while read FILE; do \
    awk -v pat="$rep" -v RS=">" '$0~pat{printf ">"; print; exit}' $FILE > ./data/annotations/rep_prot_sequences/$rep.faa ;\
  done;\
  find ./data/annotations/bakta/ -type f -iname "*.faa" -exec grep -ln "$mem" {} + | \
  while read FILE; do \
    awk -v pat="$mem" -v RS=">" '$0~pat{printf ">"; print; exit}' $FILE > ./data/annotations/mem_prot_sequences/$mem.faa ;\
  done;\
done

#RENAME all sequences to foldseek hit accession AND the contig it came from (multiple contigs may have the same protein hit, so for the sake of keeping track add in the contig name)
for i in ./data/annotations/foldseek/*_annot.txt; do \
  c=${i##*/}
  c=${c%%_besthits*}
  awk -F "\t" 'NR>1{print $1, $2}' $i | while read bk_prot acc; do \
  #move and rename header
  mv data/annotations/mem_prot_sequences/${bk_prot:0:12}.faa data/annotations/mem_prot_sequences/${acc}_$c.faa
  sed -i "s/${bk_prot:0:12}*/${acc}/g" data/annotations/mem_prot_sequences/${acc}_$c.faa ;\
  mv data/annotations/rep_prot_sequences/${bk_prot:0:12}.faa data/annotations/rep_prot_sequences/${acc}_$c.faa
  sed -i "s/${bk_prot:0:12}*/${acc}/g" data/annotations/rep_prot_sequences/${acc}_$c.faa ;\
done;\
done

##create a3m (query centered MSA) using the representatives (query) as clusters.
#for each representative in a cluster, I have to make a db of the representative and its clusters.

#GET STRUCTURES INTO RESPECTIVE CLUSTER FOLDERS
#first make folders for each cluster, and copy representative structure into its own folder
cat ./data/clustering/clu.tsv | awk '{print $1}' | while IFS="+" read x y; do echo $x; done | sort -u | while read x; do \
  mkdir ./data/clustering/clusters/${x%.pdb} ./data/clustering/clusters/${x%.pdb}/rep ./data/clustering/clusters/${x%.pdb}/members;\
  find ./data/structures/ -type f -iname "$x" -exec cp -t ./data/clustering/clusters/${x%.pdb}/rep  {} +; done
#copy member structures into respective cluster folders
cat ./data/clustering/clu.tsv | sort -u | while read x y; do \
  find ./data/structures/ -type f -iname "$y" -exec cp -t ./data/clustering/clusters/${x%.pdb}/members {} +; done 

#Use foldseek as a lookup table to determine which cluster representatives are identifiable phage related proteins (keyword phage/head/tail/capsid/baseplate)
#move phage related protein reps into their own folder
mkdir data/clustering/clusters_phage_related

for i in data/clustering/clusters/*; do \
  c=${i##*/} ;\
  c=${c:0:12} ;\
  #get files (contig fs annotations) that contain our representative, awk for col 3 containing keywords (phage/head/tail/capsid/baseplate) and get the name of the original protein (from col 1). Use those names to copy the cluster folders into a new folder (which will contain representative clusters that are annotated to be proteins that can be used for phage identification). 
  find ./data/annotations/foldseek/ -type f -iname "*_annot.txt" -exec grep -hr "$c" {} + | awk -F "\t" '{print $1, $2, $3}' |\
    while read x y z; do echo -e $x"\t"$y"\t"$z; done | \
      awk -F "\t" '$3~/phage|head|tail|capsid|baseplate/ {print $1, "\t", $2, "\t", $3}' | \
          while IFS="\t" read x y; do \
            find ./data/clustering/clusters/rep_structures -type f -iname "$x*" -exec cp {} ./data/clustering/clusters_phage_related/rep_structures/ \; ;\
            find ./data/clustering/clusters/ -maxdepth 1 -iname "$x*" -exec cp -r {} ./data/clustering/clusters_phage_related/ \; ;\
            done;\
  done

#Shorter version (gets cluster reps that contain our keyword and their structures)
for i in ./data/annotations/foldseek/*_annot.txt; do \
  awk -F "\t" 'NR>1 && $3~/phage|head|tail|capsid|baseplate/ {print $1, "\t", $2}' $i ; done | while read x y; do \
      find ./data/clustering/clusters/ -maxdepth 1 -iname "$x*" -exec cp -r {} ./data/clustering/clusters_phage_related/ \; ;\
      done

#If the rep seqeuence is in the members, remove it
for i in ./data/clustering/clusters_phage_related/*; do \
  RM=$(find $i/rep/ -type f -iname "*.pdb")
  rm $i/members/${RM##*/} ;\
done

#Get aa sequences of reps (phage related)
for i in ./data/clustering/clusters_phage_related/*; do \
  rep=${i##*/} ;\
  rep=${rep:0:12} ;\
  find ./data/annotations/bakta/ -type f -iname "*.faa" -exec grep -ln "$rep" {} + | \
  while read FILE; do \
    awk -v pat="$rep" -v RS=">" '$0~pat{printf ">"; print; exit}' $FILE > ./data/clustering/clusters_phage_related_rep_seqs/$rep.faa ;\
  done; done

#Get lookup table lines of our reps
for i in data/clustering/clusters_phage_related_rep_seqs/*.faa; do \
  x=${i##*/} ;\
  x=${x%.faa} ;\
  find ./data/annotations/foldseek/ -type f -iname "*_annot.txt" -exec grep -h "$x" {} +
done
```

Perform QC-MSA using representatives as the query, and its members as targets.

```bash
#Create fs databases
for i in data/clustering/clusters_phage_related/*; do \
  foldseek createdb $i/rep $i/querydb ;\
  foldseek createdb $i/members $i/memdb ;\
  foldseek search $i/querydb $i/memdb $i/aln $i/tmp -a ;\
  foldseek result2msa $i/querydb $i/memdb $i/aln $i/msa --msa-format-mode 6 ;\
  foldseek unpackdb $i/msa $i/msa_output --unpack-suffix .a3m --unpack-name-mode 0 ;\
  tools/reformat.pl $i/msa_output/*.a3m $i/msa_output/output.fasta;\
  done
```

Re-align our cluster representatives to its members to generate a TMscore.
Basically doing the same thing as above but the output sucks less to deal with for what I'm going to use it for

```bash
for i in data/clustering//clusters_phage_related/*; do \
  foldseek easy-search $i/members $i/querydb $i/${i##*/}_fs_results.m8 tmpFolder --alignment-type 1;\
  done
```

Get list of contigs and label samples they come from

```bash
 for i in *.fasta; do echo -ne ${i%.*} "\t"; find ../../*-*/INHERIT_phages -iname "${i##*/}" | while read x; do x=${x#*-}; echo ${x%%/*}; done; done > contigs_samples_list.txt
```

Get longest mapped reads

```bash
for i in *-*/INHERIT_phages/coverage/polished_coverage/*_sort.bam; do \
  BAM=${i##*/}; \
  echo -ne "${i%%/*}\t${BAM%%_sort*}\t"; #print sample and contig name \
  echo -ne "$(grep -v ">" ${i%%/coverage*}/${BAM%%_sort*}.fasta | tr -d "\n" | wc -c)\t"; #print size \
  samtools view $i | awk '{print length($10)}' | sort -n | tail -n 1; #print longest alignment \
done > ./phage_longest_align.txt
```
