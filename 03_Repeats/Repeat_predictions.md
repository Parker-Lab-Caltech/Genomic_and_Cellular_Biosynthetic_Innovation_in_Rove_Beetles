## **Identification and Classification of Repetitive Elements in rove beetle genomes**

Sheila Kitchen, 2/25/19

We used to methods to identify repeats in the rove beetle genomes. The first was based on the genome assembly (RepeatModeler, MITEtracker and RepeatMasker) and the second was based on the raw reads (dnaPipeTE and RepeatExplorer).

Below is how these different analyses were performed on *Dalotia*.

## Genome assembly-based predictions
### 1. [RepeatModeler v1.0.11](http://www.repeatmasker.org/RepeatModeler/) - used to identify repeat families *de novo*

```
##Make database for RepeatModeler and RepeatMasker
/central/groups/Parker_lab/tools/RepeatModeler-open-1.0.11/BuildDatabase -engine ncbi \
-name "$NAME" \
/central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/8-polish/Dcor_pilon_SPRITE.fasta

#RepeatModeler
/central/groups/Parker_lab/tools/RepeatModeler-open-1.0.11/RepeatModeler \
 -engine ncbi -pa 5 -database $NAME
```

Found 408 consensus repeat families, 314 of which were unclassified ("unknown").

I further attempted to classify the repeat families using CENSOR:
```
/central/groups/Parker_lab/tools/censor-4.2.29/bin/censor.ncbi ./${NAME}-families.fa -s -lib inv '-p 6'

#pull out the matches with the highest score for each family
sort -k 1,1 -k 11,11 -r -n ${NAME}-families.fa.map | sort -u  -k1,1 > ${NAME}_cen_uniq.map
#convert to tab-delimited
sed "s/ \+/\t/g" ${NAME}_cen_uniq.map > out
mv out ${NAME}_cen_uniq.map
```

### 2. Filter the TEs to remove "real" gene fragments
Here I filtered the consensus library following similar methodology in [Petersen *et al.* 2019](https://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-018-1324-9#Sec14). First the library was searched against a local database of beetle genomes with blastx.
```
library="${NAME}-families.fa"
blastdbcmd="/central/groups/Parker_lab/tools/ncbi/ncbi-blast-2.7.1+/bin/blastdbcmd"
blastoutfile="${NAME}_repeats_blastReport.txt"
db="/central/groups/Parker_lab/blast_databases/protDB"
fastagrep="/central/groups/Parker_lab/tools/fastagrep.pl"

echo "## Blast beetle database..."

#blast repeat families against beelte protein database
/central/groups/Parker_lab/tools/blastx -query $library -db $db -evalue 1e-5  \
-max_hsps 1 -num_threads 6 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
-out $blastoutfile

#get unique hits
echo "## Filtering for top scoring hit"
grep -v '^#' $blastoutfile \
       | sort -k 1,1n -k12,12gr -k 11,11g \
       | sort -u  -k1,1 > uniq_beetle.blast.txt
echo "$(wc -l < uniq_beetle.blast.txt) hits for all elements"

#pull out headers from the repeat model families file
grep ">" ${NAME}-families.fa > ${NAME}-head.txt

sed 's/>//' ${NAME}-head.txt > ${NAME}_head2.txt

sed 's/\s.*$//' ${NAME}_head2.txt > ${NAME}_head3.txt

rm ${NAME}-head.txt
rm ${NAME}_head2.txt

module load gcc/9.2.0
module load R/4.0.2

Rscript ../../repeatTable_merge.R ${NAME}_head3.txt uniq_beetle.blast.txt ${NAME}_cen_uniq.map
```

I merged the classifications from RepeatClassifier with the CENSOR annotation and compared them at the large superfamily level, as in a LTR transposon with RepeatClassifier was a LTR transposon with Censor. Repeat families were removed if they matched "genuine" proteins. Families that were originally classified as "unknown" and resulted in a blast hit for a known TE using a list of keywords were retained. Other predicted "unknown" TEs without a blast hit or annotation recovered using CENSOR or RepeatClassifier were removed.

See https://github.com/mptrsen/mobilome/tree/master/code for more details.

393 filtered repeat families remained.

```
fastagrep="/central/groups/Parker_lab/tools/fastagrep.pl"
library="${NAME}-families.fa"

# keep these from the repeat library, remove the rest
echo "## Filtering repeat library"
cd ./${NUM}_${SPECIES}/repeatModel
$perl $fastagrep -f /central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/repeatModel/keep-these-headers.txt \
$library > filtered-${NAME}-library.fa
echo "$(grep -c '>' filtered-${NAME}-library.fa) sequences in filtered library"
```

## 3. [Mite-Tracker](https://github.com/INTABiotechMJ/MITE-Tracker)
Identify candidate non-autonomous miniature inverted transposable elements (MITEs) in the genome.
```
python3 -m MITETracker \
-g /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/8-polish/Dcor_pilon_SPRITE.fasta \
 -w 5 -j Dcor_MT
```

MITE Tracker identified 101965 candidate MITEs in 32 clusters/families.

Classify the 32 families using RepeatClassifier (part of RepeatModler software)
```
/central/groups/Parker_lab/tools/RepeatModeler-open-1.0.11/RepeatClassifier \
-consensi /central/groups/Parker_lab/tools/MITE-Tracker/results/Dcor_MT/families_nr.fasta \
-engine ncbi
```
19 were classified and 13 remained unclassified.

### 4. Combine the custom *Dalotia* TE library with MITE Tracker and RepBase libraries.
I followed the guidelines suggested here (https://github.com/umd-byob/presentations/tree/master/2015/0324-RepeatMasker-RepeatModeler) and here (http://blaxter-lab-documentation.readthedocs.io/en/latest/repeatmasker.html).

a. Pull out the Coleoptera sequences (n= 276) in the RepDatabase (Dfam_Consensus-20181026, RepBase-20181026):

```
queryRepeatDatabase.pl -species "Coleoptera" > repeatmasker.Coleoptera.fa
```

b. Combine the Coleoptera sequences with the Dalotia custom library and MITE Tracker families, then create clusters using vsearch:

```
#combine files
cat /central/groups/Parker_lab/tools/RepeatMasker/util/repeatmasker.Coleoptera.fa \
/central/groups/Parker_lab/tools/MITE-Tracker/results/${NAME}_MT/families_nr.fasta.classified \
/central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/repeatModel/filtered-${NAME}-library.fa \
 > ${NAME}.ALL.fa

#cluster families
echo "## Cluster repeat libraries"
/central/groups/Parker_lab/tools/MITE-Tracker/vsearch-2.7.1/bin/vsearch \
--cluster_fast ${NAME}.ALL.fa --threads 4 --strand both \
--clusters ${NAME}.clust_temp --iddef 1 --id 0.8 --uc out.uc --centroids ${NAME}.repeatlib.fa
```

701 total TE families before clustering and afterwards left with 647 repeat families/clusters.

### 5. [RepeatMasker v 4.07](http://repeatmasker.org/RMDownload.html)
Use the filtered repeat library to mask the genome. In this case I used softmasking, where the repeats are lowercase and everything else is uppercase letters.
```
/central/groups/Parker_lab/tools/RepeatMasker/RepeatMasker -engine ncbi -pa 4 -s -xsmall \
-lib /central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/repeatMask/${NAME}.repeatlib.fa \
-dir ./RM_output -poly -gff \
/central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/8-polish/Dcor_pilon_SPRITE.fasta

#Build Summary for the .out file from RepeatMasker since the default assumes that it is for Homo sapiens
#create a tab-delimited file that contains the size of each chromosome/scaffold/contig in the fasta file you used
samtools faidx /central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/8-polish/Dcor_pilon_SPRITE.fasta

cp /central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/8-polish/Dcor_pilon_SPRITE.fasta.fai \
/central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/8-polish/Dcor_pilon_SPRITE.fasta.tsv

#summary script
perl /central/groups/Parker_lab/tools/RepeatMasker/util/buildSummary.pl  -species ${NAME} \
-genome /central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/8-polish/Dcor_pilon_SPRITE.fasta.tsv \
-useAbsoluteGenomeSize ./RM_output/Dcor_pilon_SPRITE.fasta.out > ./RM_output/${NAME}.out.detailed

#calculate the Kimura substitution level from the consensus sequence (history of TE accumulation in genome), CpG adjusted
perl /central/groups/Parker_lab/tools/RepeatMasker/util/calcDivergenceFromAlign.pl -s ./RM_output/${NAME}.divsum ./RM_output/Dcor_pilon_SPRITE.fasta.cat.gz
```

**Summary Table for Dalotia**

|TE category        |     sequence length|  % sequence|
|---------------------|-----------------:|------:|-------:|
|   SINEs           | 30213 | 0.00% |
|   LINEs           |  564174 | 0.13%  |
|   LTR             |       417957| 0.21% |
|   DNA             |     1794619| 0.92% |
| RC Helitron | 2951208| 1.51 %|
|   Satellites      |         65994484 | 33.78% |
|   Simple Repeats  |  1985607| 1.02% |
|   Low Complexity  |    486756| 0.25% |
|   Small RNA       |     25617  | 0.00% |
| unknown | 8796455| 4.50|
| TOTAL |83054450 | 42.52%|

## Raw read based predictions
###1. [dnaPipeTE](https://github.com/clemgoub/dnaPipeTE)

For this tool, we used a genome coverage of 25% based on the final genome assembly for each species. We acknowledge that there are difference in the assembly size to the actual genome sizes based on flow estimates that were not accounted for in this analysis. Because we applied a universal genome coverage but uncertainty exist around the genome size estimates of the rove beetles in particular, we likely assembled differing proportion of read per species and started to capture more non-repetitive sequences in the process.

```
while read name num reads size; do
        # make output directory
        mkdir ${name}_0.25_output

        #change the sample IDs, cutadapt v1.18
        cutadapt -a AGATCGGAAGAGC -m 50 -q 15 -j 0 \
        -o /central/scratchio/${user}/${name}_filtered_R1.fastq \
        /central/groups/Parker_lab/raw_genome_data/${name}/${reads}

        # run dnaPipeTE
        python3 ./dnaPipeTE.py -input /central/scratchio/${user}/${name}_filtered_R1.fastq \
        -output ./${name}_0.25_output -cpu 6 -genome_size ${size} -genome_coverage 0.25 -sample_number 2

done < samples.txt

```

To ensure that our calculations were based on presumably only repeat sequences, we filtered out non-repetitive Trinity assembled genes from the results using blast homology against gene predictions from multiple beetles and Dalotia:
```
# Dalotia
while read name reads size; do
        cd ./${name}_0.25_output

        #blast database
        /central/groups/Parker_lab/tools/blastx -query Trinity.fasta -db /central/groups/Parker_lab/final_genomes/Dalotia/Dcor_assembly_v2_200326/Dcor_assembly_v2_200326.protein.fasta -evalue 1e-5  \
        -max_hsps 1 -num_threads 4 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
        -out ${name}_Dcor_blastReport.txt
        cat  ${name}_Dcor_blastReport.txt | sort -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge >  ${name}_Dcor_blastx_results.txt
        rm ${name}_Dcor_blastReport.txt
        cd ../
done < samples.txt

# other beetle genomes
while read name reads size; do
       cd ./${name}_0.25_output

       #blast database
       /central/groups/Parker_lab/tools/blastx -query Trinity.fasta -db /central/groups/Parker_lab/blast_databases/protDB -evalue 1e-5  \
       -max_hsps 1 -num_threads 4 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle' \
       -out ${name}_blastReport.txt
       cat  ${name}_blastReport.txt | sort -k1,1 -k12,12gr -k11,11g -k3,3gr | sort -u -k1,1 --merge >  ${name}_blastx_results.txt
       rm ${name}_blastReport.txt
       cd ../
done < samples.txt
```

We combined the various output tables, parsed them for "repeat" and "non-repeat" categories using classification of Dalotia's genes and wrote these new tables to file.
```
while read name reads size; do
        echo "${name}"
        cd ./${name}_0.25_output

        # merge blast results with counts table
        Rscript ../repeatTable_merge.R reads_per_component_and_annotation ${name}_blastx_results.txt ${name}_Dcor_blastx_results.txt ../Dcor_repeats.txt Counts.txt

        cd ../

done < samples.txt
```

The repeatTable_merge.R script:
```
#!/usr/bin/env Rscript

options(warn = -1)

suppressMessages(library(dplyr))
suppressMessages(library(scales))
suppressMessages(library(tidyr))

args <- commandArgs(TRUE)
args

a<- read.table(args[1], sep=" ", quote="", header=FALSE, comment.char = "#")
b<- read.table(args[2], sep="\t", quote="", header=FALSE, comment.char = "#")
c<- read.table(args[3], sep="\t", quote="", header=FALSE, comment.char = "#")
d<- read.table(args[4], sep="\t", quote="", header=TRUE, comment.char = "#")
cnt<- read.table(args[5], sep="\t", quote="", header=FALSE, comment.char = "#")

cnt<-tail(cnt, n=1)
cnt<-cnt[,2]

print("total aligned bp")
print(cnt)

#all the blast results in one table
e<- a %>%
  left_join(
    b %>%
      select(V1, V2, V11, V13),
    by=c('V3'='V1')) %>%
  left_join(
    c %>%
      select(V1, V2, V11, V13),
    by=c('V3'='V1')) %>%
  left_join(d, by=c('V2'='Gene')) %>%
  separate(.,col=V6, into=c("Class", "Subclass"),sep= "/")

colnames(e)<-c("read_counts", "aligned_bp", "trinity_id", "RM_len", "RM_anno", "RM_class", "RM_subclass","proportion_len", "Blast_Hit", "E-value", "Hit_Description", "Dcor_Hit","Dcor_E-value","Dcor_Description", "TE")

print("no repeats table")

f<-e %>%
  filter(TE == "no repeat") %>%
  mutate(RM_class=ifelse(RM_class == "", "Unknown", RM_class)) %>%
  mutate(RM_class=ifelse(RM_class == "tRNA" | RM_class == "snRNA" | is.na(RM_class), "others" , RM_class)) %>%
  mutate(RM_class=ifelse(RM_class == "LINE?", "LINE", RM_class)) %>%
  mutate(RM_class=ifelse(RM_class == "SINE?", "SINE", RM_class))

f$RM_class<-factor(f$RM_class, levels=c("LTR", "LINE", "SINE", "DNA", "RC", "rRNA", "Low_complexity", "Satellite", "Simple_repeat", "others", "Unknown"))

f %>% group_by(RM_class) %>% summarize(total_bp=sum(aligned_bp),percent_bp=total_bp/cnt)

print("repeats table")

g<-e %>%
  filter(TE == "repeat" | is.na(TE)) %>%
  mutate(RM_class=ifelse(RM_class == "tRNA" | RM_class == "snRNA" | is.na(RM_class) , "others" , RM_class)) %>%
  mutate(RM_class=ifelse(RM_class == "", "Unknown", RM_class)) %>%
  mutate(RM_class=ifelse(RM_class == "LINE?", "LINE", RM_class)) %>%
  mutate(RM_class=ifelse(RM_class == "SINE?", "SINE", RM_class))

g$RM_class<-factor(g$RM_class, levels=c("LTR", "LINE", "SINE", "DNA", "RC", "rRNA", "Low_complexity", "Satellite", "Simple_repeat", "others", "Unknown"))

g %>% group_by(RM_class) %>% summarize(total_bp=sum(aligned_bp), percent_bp=total_bp/cnt)
g %>% summarize(non_repeat_bp=cnt-sum(aligned_bp),non_repeat_percent=non_repeat_bp/cnt)

write.table(f, "Not_repeats_table.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(g, "Repeats_table.txt", sep="\t", quote=FALSE, row.names=FALSE)
```
The final tables were visually inspected for correct repeat v. non-repeat classification in cases where no homology was found with Dalotia sequences. In those cases, matches to predicted genes in the other beetle genomes that were not annotated as a transposon or "unknown" with an e-value < 1e-10 were moved to the non-repeat table. The final, filtered repeat table was tabulated by repeat class using pivot tables in Excel.

### 2. [RepeatExplorer2](https://repeatexplorer-elixir.cerit-sc.cz/)
We turned to RepeatExplorer to try to further identify the Dalotia repeat that accounted for ~60% of the aligned reads in dnaPipeTe. This tool incorporates additional repeat databases and TAREAN pipeline, an automated identification of satellite repeats based on the topology of their cluster graphs. For RepeatExplorer2, we followed the best practices for de novo repeat identification (Procedure 1) described in Novak et al. 2020 Nature Protocols publication.

First we subsampled the reads to 2 million each for our two Dalotia individuals.
```
/central/groups/Parker_lab/tools/seqtk/seqtk sample -s 10 $dir/${NUM}_${SPECIES}/*1.fastq.gz 2000000 | gzip -c > $scratch/${NUM}_${SPECIES}/${NUM}_${SPECIES}_2M_1.fastq.gz
/central/groups/Parker_lab/tools/seqtk/seqtk sample -s 10 $dir/${NUM}_${SPECIES}/*2.fastq.gz 2000000 | gzip -c > $scratch/${NUM}_${SPECIES}/${NUM}_${SPECIES}_2M_2.fastq.gz
```
The remainder of the steps were completed using the Galaxy portal following Procedure 1. Only ~2% of the reads were used in the analysis due to the available RAM on the Galaxy portal. In both samples, ~60% of the reads belonged to a single satellite we called Dcor-Sat1. This satellite was the same repeat identified in dnaPipeTE. It was also found by RepeatModeler, making up about 33% of the predicted repeats in the genome-based repeat prediction.  

### Identify and quantify the abundant satellite Dcor-Sat1 in the Dalotia genome

To do this, we created bed files of the exons, intron and intergenic sequence in the Dalotia genome using the gene GFF file.
```
awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4,$5}' Dcor_assembly_v3_220711.combined.gff3 > Dcorv3_exon.bed
sortBed -i Dcorv3_exon.bed > Dcorv3_exon_temp.bed
mv Dcorv3_exon_temp.bed Dcorv3_exon.bed
mergeBed -i Dcorv3_exon.bed > Dcorv3_exon_merged.bed

awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4,$5}' Dcor_assembly_v3_220711.combined.gff3 > Dcorv3_gene.bed
sortBed -i Dcorv3_gene.bed > Dcorv3_gene_temp.bed
mv Dcorv3_gene_temp.bed Dcorv3_gene.bed
subtractBed -a Dcorv3_gene.bed -b Dcorv3_exon_merged.bed > Dcorv3_intron.bed

awk 'BEGIN{OFS="\t";} {print $1, "0", $2}' Dcor_assembly_v3_220711.combined.fasta.fai > Dcorv3.idx
sortBed -i Dcorv3.idx | cut -f 1,3 - > Dcorv3.idx.temp
mv Dcorv3.idx.temp Dcorv3.idx
complementBed -i Dcorv3_gene.bed -g Dcorv3.idx > Dcorv3_intergenic.bed
```

We extracted Dcor-Sat1 overlaps from the different bed files using the RepeatMasker GFF file. Ex. for introns:
```
bedtools intersect -a Dcor_assembly_v3_220711.repeats.gff -b Dcorv3_intron.bed -wo | grep "rnd-5_family-549" > repIntersect.intron
```
Then, we calculated the base pair coverage of the Dcor-Sat1 overlap to the different genomic features.
```
awk 'BEGIN{OFS="\t";}{s+=$16}END{print s}' repIntersect.intron
```

Next, we ran [RepeatProfiler](https://github.com/johnssproul/RepeatProfiler) on the Dcor-Sat1 consensus sequence against all species to get a sense of how many other species share this satellite. Satellites can be species or lineage specific and vary in abundance among those that share the repeat. In this case, the Dcor-Sat1 satellite was concatenated 4X to provide additional alignment length for the 150bp paired end reads per recommendation of the tool developers. We randomly sampled 5 million reads for each species.
```
dir=/central/groups/Parker_lab/raw_genome_data
scratch=/central/scratchio/sak3097/RP_reads

mkdir ./RepeatProfilerData

while read NUM SPECIES; do
        # subset raw reads, 5 million each
        mkdir $scratch/${NUM}_${SPECIES}
        mkdir $scratch/${NUM}_${SPECIES}/output
        /central/groups/Parker_lab/tools/seqtk/seqtk sample -s 10 $dir/${NUM}_${SPECIES}/*1.fastq.gz 5000000 | gzip -c > $scratch/${NUM}_${SPECIES}/${NUM}_${SPECIES}_5M_1.fastq.gz
        /central/groups/Parker_lab/tools/seqtk/seqtk sample -s 10 $dir/${NUM}_${SPECIES}/*2.fastq.gz 5000000 | gzip -c > $scratch/${NUM}_${SPECIES}/${NUM}_${SPECIES}_5M_2.fastq.gz

done < genome.list

# run repeatprofiler on each read set
./repeatprof profile -p ./Dcor_sat1.fasta $scratch/RP_reads -o ./RepeatProfilerData -t 4 -indel 0.5
```

Last, we ran RepeatMasker using the Dcor-Sat1 consensus sequence (4X) and 5 million read subset of Dalotia from above to estimate the Kimura's divergence. This was to investigate when it appeared and when it became abundant. 
```
#Dalotia 1
mkdir RM_output

/central/groups/Parker_lab/tools/seqtk/seqtk seq -a /central/scratchio/sak3097/RP_reads/1_Dalotia_5M_1.fastq.gz > /central/scratchio/sak3097/1_Dalotia/1_Dalotia_5M_1.fasta

/central/groups/Parker_lab/tools/RepeatMasker/RepeatMasker -engine ncbi -pa 4 -s -xsmall \
-lib Dcor_sat1.fasta \
-dir ./RM_output -poly -gff \
/central/scratchio/sak3097/1_Dalotia/1_Dalotia_5M_1.fasta

perl /central/groups/Parker_lab/tools/RepeatMasker/util/buildSummary.pl  -species Dcor_sat \
-genome /central/scratchio/sak3097/1_Dalotia/1_Dalotia_5M_1.fasta.tsv \
-useAbsoluteGenomeSize ./RM_output/1_Dalotia_5M_1.fasta.out > ./RM_output/Dcor_sat.out.detailed

#calculate the Kimura substitution level from the consensus sequence (history of TE accumulation in genome), CpG adjusted
perl /central/groups/Parker_lab/tools/RepeatMasker/util/calcDivergenceFromAlign.pl -s ./RM_output/Dcor_sat.divsum ./RM_output/1_Dalotia_5M_1.fasta.cat.gz

```
