## **Gene annotation of the Dalotia genome assembly**

Sheila Kitchen, 4/29/19

The predicted genes (evm_update2.fasta) were split into 80 chunks using fastasplit from exonerate version 2.2.0.
```
fastasplit -f evm_update2.fasta -o ./genome_split -c 80
```

These chunks were searched against NCBI nr, UniProt, pfam, merops, and cazy databases for homology using an array batch script on the Caltech HPC.

An example script:
```
#!/bin/bash
#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "nr_bl"   # job name
#SBATCH --output "nr_bl_%A-%a.out"
#SBATCH --array=0-9

user=sak3097

/central/groups/Parker_lab/tools/blastp -query ../genome_split/evm_update2.fasta_chunk_$SLURM_ARRAY_TASK_ID \
-db /central/groups/Parker_lab/tools/ncbi/nr/nr \
-evalue 1e-5 -max_target_seqs 5 -max_hsps 1 -num_threads 6 \
-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-out /central/scratchio/$user/blast_nr/Dcor_nr.$SLURM_ARRAY_TASK_ID.out
```
The domtblout for the pfam and cazy output files were further filtered using the cath-resolve-hits tool with:
```
/central/groups/Parker_lab/tools/cath-tools-0.16.2/cath-resolve-hits.ubuntu14.04 \
--input-format hmmer_domtblout \
--worst-permissible-evalue 1e-05 \
--hits-text-to-file Dcor.0.resolve Dcor_pfam.0.out
```

The blast reports were filtered for the top blast hit as follows:
```
cat  all_nr_results.txt | sort -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge >  top_nr_results.txt
```

The output files from each chunk were combined from each database separately first, and then combined into one table using an R script.
```
library(dplyr)
library(scales)
library(tidyr)

setwd('D:/Caltech/Beetle Genomes/Assemblies/1_Dalotia/gene_annotation')

#load the various results files
prot = read.table ('protein.head')
sp = read.table('top_sprot_results.txt', header=FALSE, sep = "\t",quote = "")
ncbi = read.table('top_nr_results.txt', header=FALSE, sep = "\t", quote = "")
merops= read.table('top_merops_results.txt', header=FALSE, sep = "\t", quote = "")
mh=read.table('merops.head', header=FALSE, sep = "\t", quote = "")
egg = read.table('all_eggnog_results.txt', header=FALSE, sep = "\t",quote = "")

#all the blast results in one table
b<-prot %>%
  left_join(
    sp %>%
      select(V1, V2, V11, V13),
    by='V1') %>%     
  left_join(
    ncbi %>%
      select(V1, V2, V11, V13),
    by='V1') %>%
  left_join(
    merops %>%
      select(V1, V2, V11),
    by='V1')%>%
  left_join(
    egg %>%
      select(V1, V2, V3, V5,V7,V8,V9,V10,V12,V21,V22),
    by='V1')

colnames(b)<-c("Protein_ID","Swiss-Prot_Hit","Sp_E-value","Swiss-Prot_Description",
               "NCBI_hit", "nr_E-value","NCBI_Description",
               "MEROPS_hit", "MEROPS_E-value",
               "seed_eggNOG_ortholog", "eggnog_E-value","Taxonomic_scope", "GO_terms",
               "KEGG_enzyme", "KEGG_orthology", "KEGG_map", "KEGG_reaction",
               "COG_functional_categories", "eggNOG_functional_description")
head(b)

#convert 1 column into 2, separating by the dash
mh2<-separate(mh,V1,c("MEROPS_hit", "MEROPS_Description"), " - ")

#merge with the merops description table
c<-b %>%
  left_join(
    mh2 %>% # join with the second file (only the first and third column)
      select('MEROPS_hit','MEROPS_Description'),
    by='MEROPS_hit')

head(c)

#change order of columns
c = c %>% select("Protein_ID","Swiss-Prot_Hit","Sp_E-value","Swiss-Prot_Description",
                 "NCBI_hit", "nr_E-value","NCBI_Description",
                 "MEROPS_hit", "MEROPS_E-value","MEROPS_Description",
                 "seed_eggNOG_ortholog", "eggnog_E-value","Taxonomic_scope", "GO_terms",
                 "KEGG_enzyme", "KEGG_orthology", "KEGG_map", "KEGG_reaction",
                 "COG_functional_categories", "eggNOG_functional_description")

write.table(c, "Dcor_GeneAnnotation_combined.txt", sep="\t", quote=FALSE, row.names=FALSE)

#load in pfam and cazy results
cazy= read.table('top_cazy.txt', header=FALSE, sep = " ", quote = "")
dbca=read.table('dbCAN-fam-HMMs.txt', header=FALSE, stringsAsFactors=TRUE,
                sep = "\t", quote ="")
pfam=read.table('top_pfam.out', header=FALSE, sep = " ", quote = "")
dbpf=read.table('Pfam-A.clans.tsv', header=FALSE,
                sep = "\t", quote ="")

head(cazy)
head(dbca)
head(pfam)
head(dbpf)

cazy$V2<-gsub("len=", "len_", cazy$V2)
pfam$V1<-gsub("len=", "len_", pfam$V1)

#cazy table- reorder with gene id first
d<-cazy %>%
  select(V2,V1,V3,V4,V5,V6,V7)

colnames(d)<-c("V1","V2","V3","V4","V5","V6","V7")

#merge with prot table
e<-prot %>%
  full_join(
    d %>%
      select(V1, V2, V4, V6),
    by='V1')

#change header of the cazy description table
colnames(dbca)<-c("V2","Description","Vx")

#merge with cazy descriptions
f<-e %>%
  left_join(
    dbca %>%
      select('V2','Description'),
    by='V2')

head(f)

#write to file
write.table(f, "Dcor_GeneAnnotation_Cazy.txt", sep="\t", quote=FALSE, row.names=FALSE)

#pfam table
g<-pfam %>%
  select(V1,V2,V4,V6)

colnames(g)<-c("V1","V2","V3","V4")

h<-prot %>%
  full_join(
    g %>%
      select(V1, V2, V3, V4),
    by='V1')
head(h)

colnames(dbpf)<-c("Pfam_id","Clan_id","Clan_name","V2", "Domain_Description")
head(dbpf)

i<-h %>%
  left_join(
    dbpf %>%
      select('V2',"Pfam_id","Clan_id","Clan_name",'Domain_Description'),
    by='V2')

head(i)

colnames(i)<-c("gene_id","Domain_id","sequence_range","e-value",
               "Pfam_id","Clan_id","Clan_name","Domain_Description")

#write to fiie
write.table(i, "Dcor_GeneAnnotation_pfam.txt", sep="\t", quote=FALSE, row.names=FALSE)
```
Gene annotation was decided in order where first it takes on the Swiss-Prot annotation if evalue < 1e-10, otherwise the NCBI annotation is taken if evalue < 1e-10. If the NCBI evalue is larger than 1e-10, then the gene takes on the *T. castaneum* annotation from the eggnog database. If Tcas annotation is not available then it becomes a "hypothetical protein".

Make a tab-delimited annotation file called dcor_annotations.txt as follows with no spaces in the description:

```
{mRNA name} {protein} {annotation_here}
```

Add annotation to the gff file using [GAG](https://genomeannotation.github.io/GAG/)- must be run with python 2.7.15!!:
```
#replace "pilon" and "=" in the gff
sed 's/_pilon_pilon//g' evm_update2.gff3 > out.gff3

sed 's/len=/len_/g' out.gff3 > out2.gff3

mv out2.gff3 evm_update_pasa2x.gff3

#replace "pilon" and "=" in the genome file
sed 's/_pilon_pilon//g' Dcor.fna.masked_2.fa > out.fa

mv out.fa Dcor.fna.masked_2.fa

#remove NCBI species from which the annotation comes from
sed 's/\[[^]]*\]$//' test > test2

#remove trailing _ in the dcor_annotations.txt
sed 's/_*$//g' test2 > test3

mv test3 dcor_annotations.txt

#run GAG
python gag.py --fasta /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/Dcor.fna.masked_2.fa \
--gff /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/evm_update_pasa2x.gff3 \
--fix_terminal_ns \
--fix_start_stop \
--anno dcor_annotations.txt \
--out Dcor_gag_output

#combine the ignored gff back with the genome.gff3
cat genome.gff genome.ignored.gff > Dcor.gff3
```

Add annotation to the gene models using:
```
#remove everything but name of protein/mRNA
sed 's/\ ID=evm.model.*//g' genome.proteins.fasta > test2
sed 's/protein|/Dcor_/' test2 > test3
mv test3 Dcor.protein.fasta

sed 's/\ID=evm.model.*//g' genome.mrna.fasta > test4
sed 's/>evm/>Dcor_evm/' test4 > test5
mv test5 Dcor.mrna.fasta

#add annotations
cut -f 1,3 dcor_annotations.txt > Names_mrna.txt

perl addAnno.pl Dcor.mrna.fasta Names_mrna.txt > Dcor_assembly_v1_190525.mrna.fa

perl addAnno.pl Dcor.protein.fasta Names_mrna.txt > Dcor_assembly_v1_190525.protein.fa
```
