# **Illumina Genome Assembly Guide**

Sheila Kitchen, 5/2/2019

## Export tools to you path
export PATH="/central/groups/Parker_lab/tools/anaconda/bin:$PATH"

## 1.	Read Quality Assessment
#### [FastQC v 0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
FastQC is used to check the quality of the raw sequencing data prior to read filtering.
```
fastqc -o ./fastqc_report -f fastq -t 4 \
/central/groups/Parker_lab/raw_genome_data/31_*/*_1.fastq \
/central/groups/Parker_lab/raw_genome_data/31_*/*_2.fastq
```
-o = output directory
-f = fastq reads (can be fastq, bam or sam)
-t = threads

#### [cutadapt v1.18](https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage)
Now trim the Illmina adapters, remove bases that have a quality score below 15 from the 3' end, and those resulting reads that are now shorter than 50 bp.
```
#load correct python module
module load python3/3.7.0

#make output directories in scratch
mkdir /central/scratchio/sak3097/beetle/31/
mkdir /central/scratchio/sak3097/beetle/31/filtered

#run cutadapt on the raw reads and write filtered reads to scratch
/central/groups/Parker_lab/tools/cutadapt/bin/cutadapt \
-a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 50 -q 15 -j 0 \
-o /central/scratchio/sak3097/beetle/31/filtered/31_filtered_R1.fastq \
-p /central/scratchio/sak3097/beetle/31/filtered/31_filtered_R2.fastq \
/central/groups/Parker_lab/raw_genome_data/31_*/*_1.fastq \
/central/groups/Parker_lab/raw_genome_data/31_*/*_2.fastq
```
-a, -A = adapter
-o = output
-p = paired output
-m = minimum length
-q = quality-cutoff before adapter removal
-j 0 = automatically detect number of cores available to use

## 4. Short-Read Illumina Assembly
Use the adapter-trimmed and filter reads from above (cutadapt).

#### Contig Assembly- [MEGAHIT v1.1.3](https://github.com/voutcn/megahit)
MEGAHIT is metagenomic assembler of large and/or complex sequencing reads using de Bruijn graph.

```
module load python/2.7.15

#contig assembly
megahit --k-list 21,29,39,59,79,99,119 -m 0.5 -t 10 \
-o mh_run1 -1 /central/scratchio/sak3097/beetle/31/filtered/31_filtered_R1.fastq\
 -2 /central/scratchio/sak3097/beetle/31/filtered/31_filtered_R2.fastq
```
--k-list = comma separated list of kmer sizes  
-m = max memory (fraction of the total memory)  
-t = threads  
-o = output directory  
-1/-2 = comma-separated list of fasta/q paired-end #1/#2 files, can be compressed gz/bz2  

To get a quick summary of the assembly statistics, use the utility script from Abyss.

```
abyss-fac -t 0 final.contigs.fa
```
-t = minimum contig/scaffold length

#### Contamination Screening of the contigs -
##### [Blobtools v1.0](https://blobtools.readme.io/docs)

Blobtools is used to classify and visualize the contigs based on G/C content and read coverage.
```
#setup pointers to common files/directories
NUM=31
fastq=/central/scratchio/sak3097/beetle/1/filtered/?.fastq
assembly=/central/groups/Parker_lab/SheilaK/genome_assemblies/31_Sceptobius_lat/megahit/mh_run1/final.contigs.fa
dir=/central/scratchio/sak3097/beetle/$NUM/mbb
nt=/central/groups/Parker_lab/tools/ncbi/nt/nt
out=/central/groups/Parker_lab/SheilaK/genome_assemblies/31_Sceptobius_lat/megahit/filtered_seqs

#map reads to genome assembly
bwa index $assembly
bwa mem -t 6 $assembly $fastq | samtools view -b - > /central/scratchio/sak3097/beetle/$NUM/mbb/31_megahit.bam

#load conda environment with blobtools installed
source activate python2

#calculate read coverage per contig
blobtools map2cov -i $assembly -b /central/scratchio/sak3097/beetle/$NUM/mbb/31_megahit.bam

#blast contigs to NCBI nt
export BLASTDB=$BLASTDB:/central/groups/Parker_lab/tools/ncbi/nt/
blastn -query $assembly -db $nt -evalue 1e-25 \
-outfmt '6 qseqid sseqid sscinames staxids bitscore' \
-out $dir/31.cul10.out -culling_limit 10 -max_hsps 1 -num_threads 6

cut -f 1,4,5 $dir/31.cul10.out > $dir/31.cul10_cut.out

#create blobtools database
blobtools create -i $assembly -c 31_megahit.bam.cov -t $dir/31.cul10_cut.out -o $dir/blob_31

#create table of taxonomy
blobtools view -i $dir/blob_31* -r superkingdom -r phylum -r order --hits

#plot the results
blobtools plot -i $dir/blob_31* -x bestsum -r phylum --format png -o 31_blobplot
```
The resulting plot of the contigs based on phylum classification. The low G/C content, high coverage contigs are repetitive content.

Next, pull out the contigs that are classified as bacterial based on the blast report. This is a **conservative** approach where only those contigs that matched with high homology to the NCBI nr database are removed.
```
#create list of contigs matching bacteria
awk '$6=="Bacteria"' blob_31.blobDB.table.txt | cut -f1 > $out/bacteria.contigIDs.txt

#pull out contigs matching bacteria (n= 1,503)
blobtools seqfilter -i $assembly -l $out/bacteria.contigIDs.txt -o $out/${NUM}_bacteria
```

To create filtered scaffolds without the bacterial contigs:
```
#keep all scaffolds that are not bacteria (n= 99,494)
blobtools seqfilter -i $assembly -l $out/bacteria.contigIDs.txt \
-v -o $out/${NUM}_contigs
```
#### Check for redundancy with [Redundans v0.14a](https://github.com/lpryszcz/redundans), complied Oct 29, 2018
```
module load python/2.7.15

NUM=37
fastq=/central/scratchio/sak3097/beetle/${NUM}/filtered/${NUM}_filtered_R?.fastq

/central/groups/Parker_lab/tools/redundans/redundans.py -v -i $fastq \
-f ./bacterial_filter/${NUM}_contigs.final.contigs.filtered.fna -o ./red_run_${NUM} --iters 4 -t 10 -m 100

/central/groups/Parker_lab/tools/redundans/redundans.py -v -i $fastq -f ./bacterial_filter/${NUM}_contigs.final.contigs.filtered.fna -o ./red_run_${NUM}_dalRef \
--limit 0.5 --reference /central/groups/Parker_lab/final_genomes/Dalotia/Dcor_assembly_v1_190722.pruned_scaffolds.fa --iters 3 -t 10 -m 100

```
Remove contigs smaller than 1kbp:

```
/central/groups/Parker_lab/tools/seqtk/seqtk seq -L 1000 ./red_run_${NUM}/scaffolds.reduced.fa > ./red_run_${NUM}/scaffolds.reduced_1k.fa
```

#### Gap filling with [GapCloser](https://sourceforge.net/projects/soapdenovo2/files/GapCloser/)

This is part of the SOAPdenovo2 package and uses the same config file format.
```
/central/groups/Parker_lab/tools/redundans/bin/GapCloser \
-a ${NUM}_sspace_blob.final.scaffolds.fasta \
-b ./config -o ${NUM}_filled.fa -l 151 -t 10
```
