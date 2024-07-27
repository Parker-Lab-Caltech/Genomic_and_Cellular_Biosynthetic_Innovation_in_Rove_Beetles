## **Gene prediction of the Dalotia genome assembly**

Sheila Kitchen, 3/13/19

### 1. *Ab inito* gene prediction

#### A. GeneMark-ES v.4.33

```
# run GeneMark-ES
#load the correct conda environment and Perl library
source activate funannotate
export PERL5LIB=/central/groups/Parker_lab/tools/anaconda/lib/site_perl/5.26.2

/central/groups/Parker_lab/tools/gm_et_linux_64/gmes_petap/gmes_petap.pl --ES --cores 4 --soft_mask \
--max_intron 300000 \
--sequence /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked

#Converted gtf output to gff3 using:
/central/groups/Parker_lab/tools/funannotate-1.5.1/util/genemark_gtf2gff3.pl genemark.gtf > genemark-ES.gff3

#Validate it will work for EVM:
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/gff3_gene_prediction_file_validator.pl genemark-ES.gff3
```

### 2. Evidence based gene prediction
#### A. BRAKER - uses RNAseq Reads
##### i. Map reads to genome with STAR

Make index of genome
```
source activate base

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/ \
--genomeFastaFiles /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked
```
Map reads to the genome (same RNAseq reads use to assembly the transcriptome below). I ended up splitting this up because of memory issues in writing out the .bam file.
```
DIR=/central/groups/Parker_lab/transcriptome_data/Dalotia/
STAR --runThreadN 10 \
--genomeDir /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/ \
--outSAMtype BAM Unsorted \
--twopassMode Basic \
--outFileNamePrefix /central/scratchio/sak3097/beetle/STAR/STAR_mapping \
--readFilesCommand zcat \
--alignIntronMax 300000 \
--readFilesIn $DIR/larvae/Dal-Larvae-RNA-9-18-18-2_S26_L003_R1_001.fastq.gz,$DIR/larvae/Dal-Larvae-RNA-9-18-18-1_S15_L002_R1_001.fastq.gz,$DIR/pupae/Dal-Pupae-RNA-9-18-18-2_S25_L003_R1_001.fastq.gz,$DIR/pupae/Dal-Pupae-RNA-9-18-18-1_S14_L002_R1_001.fastq.gz,$DIR/MA/Dal-MA-RNA-9-4-18-1_S12_L002_R1_001.fastq.gz,$DIR/MA/Dal-MA-RNA-9-4-18-2_S23_L003_R1_001.fastq.gz,$DIR/FA/Dal-FA-RNA-9-20-18-1_S11_L002_R1_001.fastq.gz,$DIR/FA//data/Dal-FA-RNA-9-20-18-2_S22_L003_R1_001.fastq.gz,$DIR/MW/Dal-MW-RNA-9-20-18-1_S13_L002_R1_001.fastq.gz,$DIR/MW/Dal-MW-RNA-9-20-18-2_S24_L003_R1_001.fastq.gz,$DIR/gland_control/gland/gland_1/gland_1_1_1.fastq.gz,$DIR/gland_control/gland/gland_1/gland_1_2_1.fastq.gz,$DIR/gland_control/gland/gland_2/gland_2_1_1.fastq.gz,$DIR/gland_control/gland/gland_2/gland_2_2_1.fastq.gz,$DIR/gland_control/gland/gland_3/gland_3_1_1.fastq.gz,$DIR/gland_control/gland/gland_3/gland_3_2_1.fastq.gz,$DIR/gland_control/control/control_1/control_1_1_1.fastq.gz,$DIR/gland_control/control/control_1/control_1_2_1.fastq.gz,$DIR/gland_control/control/control_2/control_2_1_1.fastq.gz,$DIR/gland_control/control/control_2/control_2_2_1.fastq.gz,$DIR/gland_control/control/control_3/control_3_1_1.fastq.gz,$DIR/gland_control/control/control_3/control_3_2_1.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-01_S122_L008_R1_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-01_S42_L007_R1_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-05_S126_L008_R1_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-05_S46_L007_R1_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-09_S130_L008_R1_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-09_S50_L007_R1_001.fastq.gz $DIR/larvae/Dal-Larvae-RNA-9-18-18-2_S26_L003_R2_001.fastq.gz,$DIR/larvae/Dal-Larvae-RNA-9-18-18-1_S15_L002_R2_001.fastq.gz,$DIR/pupae/Dal-Pupae-RNA-9-18-18-2_S25_L003_R2_001.fastq.gz,$DIR/pupae/Dal-Pupae-RNA-9-18-18-1_S14_L002_R2_001.fastq.gz,$DIR/MA/Dal-MA-RNA-9-4-18-1_S12_L002_R2_001.fastq.gz,$DIR/MA/Dal-MA-RNA-9-4-18-2_S23_L003_R2_001.fastq.gz,$DIR/FA/Dal-FA-RNA-9-20-18-1_S11_L002_R2_001.fastq.gz,$DIR/FA/Dal-FA-RNA-9-20-18-2_S22_L003_R2_001.fastq.gz,$DIR/MW/Dal-MW-RNA-9-20-18-1_S13_L002_R2_001.fastq.gz,$DIR/MW/Dal-MW-RNA-9-20-18-2_S24_L003_R2_001.fastq.gz,$DIR/gland_control/gland/gland_1/gland_1_1_2.fastq.gz,$DIR/gland_control/gland/gland_1/gland_1_2_2.fastq.gz,$DIR/gland_control/gland/gland_2/gland_2_1_2.fastq.gz,$DIR/gland_control/gland/gland_2/gland_2_2_2.fastq.gz,$DIR/gland_control/gland/gland_3/gland_3_1_2.fastq.gz,$DIR/gland_control/gland/gland_3/gland_3_2_2.fastq.gz,$DIR/gland_control/control/control_1/control_1_1_2.fastq.gz,$DIR/gland_control/control/control_1/control_1_2_2.fastq.gz,$DIR/gland_control/control/control_2/control_2_1_2.fastq.gz,$DIR/gland_control/control/control_2/control_2_2_2.fastq.gz,$DIR/gland_control/control/control_3/control_3_1_2.fastq.gz,$DIR/gland_control/control/control_3/control_3_2_2.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-01_S122_L008_R2_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-01_S42_L007_R2_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-05_S126_L008_R2_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-05_S46_L007_R2_001.fastq.gz,$DIR/brain/female_head_group17104FL-18-01-09_S130_L008_R2_001.fastq.gz,$DIR/brain/female_head_group/17104FL-18-01-09_S50_L007_R2_001.fastq.gz

#sort the bam files
/central/groups/Parker_lab/tools/samtools-1.9/samtools sort \
-m 7G -o /central/scratchio/sak3097/beetle/STAR_mapping_1Aligned_sorted.bam \
-T /central/scratchio/sak3097/beetle/STAR_mapping_1Aligned_temp \
--threads 10 /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/BRAKER/STAR_mapping_1Aligned.bam

/central/groups/Parker_lab/tools/samtools-1.9/samtools sort \
-m 7G -o /central/scratchio/sak3097/beetle/STAR_mapping_2Aligned_sorted.bam \
-T /central/scratchio/sak3097/beetle/STAR_mapping_2Aligned_temp \
--threads 10 /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/BRAKER/STAR_mapping_2Aligned.bam
```

##### ii. Run BRAKER v 2.1.2
Run Braker with the .bam file(s) from STAR
```
source activate braker2

braker.pl \
--workingdir=/central/scratchio/sak3097/beetle –-softmasking --species=Dalotia --cores 10 \
--genome=/central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/BRAKER/genome.fa \
--AUGUSTUS_CONFIG_PATH=/central/groups/Parker_lab/tools/anaconda/envs/braker2/config \
--AUGUSTUS_BIN_PATH=/central/groups/Parker_lab/tools/anaconda/envs/braker2/bin \
--GENEMARK_PATH=/central/groups/Parker_lab/tools/gm_et_linux_64/gmes_petap \
--BAMTOOLS_PATH=/central/groups/Parker_lab/tools/anaconda/envs/braker2/bin \
--SAMTOOLS_PATH=/central/groups/Parker_lab/tools/samtools-1.9 \
--BLAST_PATH=/central/groups/Parker_lab/tools/ncbi/ncbi-blast-2.7.1+/bin \
--PYTHON3_PATH=/central/groups/Parker_lab/tools/anaconda/envs/braker2/bin \
--bam=/central/scratchio/sak3097/beetle/STAR_mapping_1Aligned_sorted.bam,/central/scratchio/sak3097/beetle/STAR_mapping_2Aligned_sorted.bam
```

##### iii. Pulled out the high-quality gene models (>90% read coverage on the exon boundaries) from the BRAKER gff3 file using grep

```
#Lines with mRNA from file:
grep -P "transcript\t" augustus.hints.gff > mrna.txt

#Lines with statement 'CDS introns':
grep "CDS introns" augustus.hints.gff | cut -f 4 -d " " - > introns.txt

paste mrna.txt introns.txt > join.txt

sed 's/\//\t/' join.txt | awk -v OFS='\t' '{if( $11 != 0) $12 = ($10 / $11)*100}1' | sort -k12,12nr | awk -F"\t" '$12>90' |cut -f 9 - > hiq_genes.txt
```

Opened up files in Excel and matched the transcript id to the CDS intron results. Calculate the coverage of that model by dividing those introns with hints by total introns (for example 2/4= 50%). Filter the list by >=90% support by the RNAseq reads.

Created a list of the top hits and pulled them out of the gff3 file using grep with a file (-f). New gff3 file was put in as the HiQ gene models for EVM.
```
/central/groups/Parker_lab/tools/augustus-3.2.3/scripts/gtf2gff.pl < augustus.hints.gff --out=out.gff --gff3

grep -f hiq_genes.txt ./out.gff > hiq.gff

/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl

assembler-sample_mydb_pasa_holobus1.sqlite

/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/gff3_gene_prediction_file_validator.pl ./input_gffs/${NUM}_augustus.gff3
```
Convert the two Augustus outputs to EVM format:
```
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl \
augustus.hints.gff > out_aug
mv out_aug augustus.hints.gff

/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl \
hiq_genes.gff > hq_aug
mv hq_aug hiq_genes.gff
```

#### B. PASA v2.3.3- uses transcriptome assembly
The reference-guided transcriptome was assembled using Trinity (same RNAseq reads as above). The transcriptome has 128,720 transcripts. This was reduced to the longest orf for each group of transcript isoforms using TRANSDECODER in the PASA script.

```
#Load conda environment
source activate funannotate

PASAHOME=/central/groups/Parker_lab/tools/PASApipeline-v2.3.3

#clean up the transcriptome files
$PASAHOME/bin/seqclean  20190313_dalotia_many_reads_genome_guided.fasta

# PASA
$PASAHOME
/central/groups/Parker_lab/tools/PASApipeline-v2.3.3/Launch_PASA_pipeline.pl \
     -c alignAssembly.config -C -R  --MAX_INTRON_LENGTH 300000 \
     -g  /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked \
     -t 20190313_dalotia_many_reads_genome_guided.fasta.clean -T \
     -u 20190313_dalotia_many_reads_genome_guided.fasta \
     --ALIGNERS gmap,blat --CPU 8 --TRANSDECODER
```

#### C. Exonerate- protein evidence  

The Insecta protein sequences were downloaded from Uniprot and used as evidence for protein matches with at least a percent identity of 80%.
```
source activate funannotate

#split the genome up by scaffold into individual fasta files
fastaexplode -f /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked_2.fa \
-d ./

#search each chunk against the Insecta sequences in swissprot/trembl
exonerate --model p2g --showvulgar no --showalignment no --showquerygff no \
-Q protein -T dna \
--showtargetgff yes --percent 80 \
--ryo "AveragePercentIdentity: %pi\n" \
-q /central/groups/Parker_lab/tools/UniProt_db/uniprot-insecta.fasta \
-t /central/scratchio/sak3097/beetle/Dcor.fna.masked_chunk_0000000 \
--softmasktarget yes >> Dal_exo_insecta_1.out

#combine the output for all individual scaffolds
cat *.gff3 > exonerate_raw.gff3

#count the number of total hits
grep -c "START OF GFF DUMP" exonerate_raw.gff3
#12890

#convert the exonerate gff format to the evidence modeler format
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/misc/Exonerate_to_evm_gff3.pl \
exonerate_raw.gff3 > exonerate_evm.gff3

#count the number of unique hits to trembl seqs
grep "tr|" exonerate_evm.gff3 | sort -u -k1,1 | wc -l
#41

#count the number of unique hits to swissprot seqs
grep -c "sp|" exonerate_evm.gff3 sort -u -k1,1 | wc -l
#31
```

### 3. Combine gene predictions with EVidenceModeler
###### To run EVM:

Each GFF files needs a unique id in column two that matches the weights.txt.
```
#remove pipes
#change name of the 2nd column
sed 's/|/_/g' genemark-ES.gff3 > genemark

sed 's/assembler-sample_mydb_pasa_dalotia5.sqlite/pasa/' \
sample_mydb_pasa_dalotia5.sqlite.pasa_assemblies.gff3 > pasa.gff3

sed 's/|/_/g' pasa.gff3 > pasa

sed 's/AUGUSTUS/HiQ/' hiq.gff > hiq_genes.gff
```

Concatenate all the gff files into one, gene_predictions.gff3.
```
cat augustus.hints.gff hiq_genes.gff genemark-ES.gff3 > gene_prediction.gff3
```
weights.txt file:
PROTEIN   exonerate  1
TRANSCRIPT      pasa   10
OTHER_PREDICTION        HiQ     4
ABINITIO_PREDICTION      Augustus       1
ABINITIO_PREDICTION     GeneMark.hmm        1

1. Partition the genome (runMe.sh script pt 1)
```
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/partition_EVM_inputs.pl \
--genome /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked_2.fa \
--gene_predictions ./input_gffs/gene_prediction.gff3 \
--protein_alignments ./input_gffs/exonerate_evm.gff3 \
--transcript_alignments ./input_gffs/pasa.gff3 \
--segmentSize 100000 --overlapSize 10000 --partition_listing Dcorpartitions_list.out
```

2. Generate list of commands (runMe.sh script pt 1)
```
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/write_EVM_commands.pl \
--genome /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked_2.fa \
--weights /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/EVM/weights.txt \
--gene_predictions ./input_gffs/gene_prediction.gff3 \
--protein_alignments ./input_gffs/exonerate_evm.gff3 \
--transcript_alignments ./input_gffs/pasa.gff3 \
--output_file_name Dcor_evm.out \
--partitions Dcorpartitions_list.out >  Dcor_commands.list
```                                                            

3. Run commands using parallel
```
module load parallel/20180222
parallel --jobs 10 < Dcor_commands.list
```

4. Combine partitions
```
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/recombine_EVM_partial_outputs.pl \
--partitions Dcorpartitions_list.out --output_file_name Dcor_evm.out
```

5. Convert to GFF3 format
```
/central/groups/Parker_lab/tools/evidencemodeler/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
--partitions Dcorpartitions_list.out --output Dcor_evm.out \
--genome /central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeatMask/DcorRMaskOUT/Dcor.fna.masked_2.fa

find . -regex ".*Dcor_evm.out.gff3" -exec cat {} \; > Dcor_EVM.all.gff3
```


To summarize which gene predictor contributed to each gene model run summarize.py:
```
#!/usr/bin/python
import sys, os

ds = os.listdir(sys.argv[1])

for d in ds:
    fPath = sys.argv[1] + '/' + d + '/Dcor_evm.out'
    size = os.path.getsize(fPath)
    if size > 0:
        blocks = open(fPath).read().strip().split('#')[1:]
        for block in blocks:
            coords = []
            evidence = []
            for line in block.strip().split('\n')[1:]:
                if line.strip() != '' and line[0] != '!':
                    meta = line.strip().split('\t')
                    coords.append(int(meta[0]))
                    coords.append(int(meta[1]))
                    coords.sort()
                    evidence.extend([tuple(x[1:-1].split(';')) for x in meta[-1].split(',')])

            evidence = set(evidence)
            sources = set([x[1] for x in evidence])

            print d + '\t' + str(coords[0]) + '\t' + str(coords[-1]) + '\t' + ','.join([x[0] for x in evidence]) + '\t' + ','.join(sources)
```
Copy all the evm output into a new directory and run as:

```
./summarize.py ./evm_output > summary.txt
```
