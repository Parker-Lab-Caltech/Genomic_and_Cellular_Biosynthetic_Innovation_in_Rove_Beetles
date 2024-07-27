## **SPRITE assembly of the Dalotia genome**

Sheila Kitchen, 1/8/2020

### Sample preparation
92 males were used for the final sample prep. The methods can be found at: https://www.lncrna-test.caltech.edu/protocols/SPRITE_Protocol_DNA_January_2018.pdf

1. The beetles were macerated with a glass dounce in 8ml of 2mM DSG solution at room temperature. Pipette up and down to get a single cell suspension. Rock gently for 45 mins. \*Due to the volume of the starting material, we doubled the disuccinimidyl glutarate solution (DSG; made fresh by adding 16 ul of 0.5M DSG in 8 ml PBS) compared to their normal starting volumes.
2. Premix 3% paraformaldehyde solution in PBS just before use. Combine 1500 ul of 16% FA with 6.5 PBS.
3. The cells + cuticle was pelleted for 8 mins at 2500xg at room temp. Remove the DSG solution.
4. The pellet was resuspended in 8ml of PBS, then re-pelleted as before.
5. Add 8ml of 3% FA solution and pipette up and down several times. Rock gently at room temp for exactly 10 mins.
6. Immediately add 200ul of 2.5M glycine stop solution per 1ml of FA solution directly to quench the crosslink reaction. Mix well.
7. Rock gently at room temp for 5 mins.
8. Spin down cells at 4C for 4 mins at 2500xg. Discard FA solution and keep cells at 4C moving forward.
9. Resuspend pellet in 8ml of cold 1x PBS + 0.5% BSA, rock for 1-2 mins.
10. Spin down cells at 4C for 4 mins at 2500xg Discard waste.
11. Repeat wash step.
12. Resuspend in 2 ml of 1x PBS + 0.5% BSA.
13. Cells were aliquoted into microcentrifuge and spin at 4C for 5min at 2000xg. Remove supernatant without distrupting the pellet.
14. Flash	freeze in	liquid	nitrogen	and	store	pellet	at	-80C.

### Library preparation
Sofi and Elizabeth finished the library prep steps to get 8 final libraries. One of these was run on the MiSeq to assess the quality of the libraries before proceeding.

The	following	metrics were used	to	evaluate	whether	SPRITE	tagging was	successful on	the	post-sequencing	library:
1. Percentage	of	reads	with	all	tags	ligated (get_ligation_efficiency.py)

|Number of Tags| Dalotia| Human\*|
|----------|-----:|---:|
|0|0.7%|0.2%|
|1|1.3%|2.7%|
|2|4.8%|8.8%|
|3|5.8%|9.8%|
|4|9.9%|12.2%|
|5|77.5%|66.3%|

This was surprising good and better than the reported values for mice and human\* samples.

2. Percentage	 of	 chromatin that	is	interacting	with other chromatin:

This ensures that the sonication is sufficient to create different clusters sizes and does not result in the majority being non-interacting sequences.

Number of reads for each barcode (interacting molecules):

![alt text](./cluster_sizes.png "Cluster Distribution")

3. FastQC	to	QC	quality	of	reads	on	sequencer.

These different metrics gave us the confidence to proceed with the HiSeq sequencing run.

### Data Processing

#### 1. Trim the Illumina adaptors with cutadapt

```
#load correct python module
module load python3/3.7.0

#trim reads
for NUM in $(cat sample_list.txt); do
        #mkdir /central/scratchio/sak3097/beetle/SPRITE
        #mkdir /central/scratchio/sak3097/beetle/SPRITE/${NUM}
        #mkdir /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered

        cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 50 -q 15 -j 0 \
        -o /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_filtered_R1.fastq \
        -p /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_filtered_R2.fastq \
        /central/groups/Parker_lab/SheilaK/sprite_dalotia/workup3/fastqs/FT-SE5482/${NUM}_R1_001.fastq.gz \
        /central/groups/Parker_lab/SheilaK/sprite_dalotia/workup3/fastqs/FT-SE5482/${NUM}_R2_001.fastq.gz
done
```
After trimming 414 million reads remained split over 20 different sequencing libraries.

#### 2. Identify barcodes in the filtered sequences
Use the java script provided by the Guttman lab to identify the barcodes in  each library separately.
```
source activate funannotate

java -jar /groups/guttman/software/sprite-pipeline/java/BarcodeIdentification_v1.2.0.jar \
--input1 /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_filtered_R1.fastq \
--input2 /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_filtered_R2.fastq \
--output1 /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_R1.barcoded.fastq.gz \
--output2 /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_R2.barcoded.fastq.gz \
--config ./workup3/6_config.txt &> ./workup3/fastqs/${NUM}.id.log
```

Calculate the ligation efficiency.
```
python /groups/guttman/software/sprite-pipeline/python/get_ligation_efficiency.py \
/central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_R1.barcoded.fastq.gz > ./workup3/${NUM}.ligation_efficiency_6_config.txt
```
Average ligation efficiency was 75.04% with 311 million reads remaining.

Filter out the short barcodes from the first read that contains the genomic DNA.
```
python /groups/guttman/software/sprite-pipeline/python/get_full_barcodes.py \
--r1 /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_R1.barcoded.fastq.gz
```

#### 3. Align reads with full barcodes to the Dalotia genome

```
module load samtools/1.8
module load bowtie/2.3.4.1

for NUM in $(cat sample_list.txt); do
        bowtie2 -p 10 -t --phred33 -x /central/groups/Parker_lab/final_genomes/Dalotia/dalotia \
        -U /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}_R1.barcoded_full.fastq.gz \
        | samtools view -bq 20 -F 4 -F 256 - > \
        /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}.bowtie2.mapq20.bam
done
```
The average alignment rate was 67.36% of the filtered reads.

#### 4. Filter for overlap with repetitive DNA

```
source activate funannotate

for NUM in $(cat sample_list.txt); do
        /groups/guttman/software/bedtools2/bin/bedtools intersect -v \
        -a /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}.bowtie2.mapq20.bam \
        -b /central/groups/Parker_lab/final_genomes/Dalotia/Dcor.fna.repeats_pruned.gff > \
        /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}.bowtie2.mapq20_masked.bam
done
```
The reads were reduced to only 22.8% on average per library after quality and repeat filtering.
#### 5. Cluster generation
```
module load R/3.6.1

for NUM in $(cat sample_list.txt); do
        python /groups/guttman/software/sprite-pipeline/python/get_clusters.py \
        --input /central/scratchio/sak3097/beetle/SPRITE/${NUM}/filtered/${NUM}.bowtie2.mapq20_masked.bam --output ./workup3/clusters/${NUM}.clusters --num_tags 5

        Rscript /groups/guttman/software/sprite-pipeline/r/get_cluster_size_distribution.r ./workup3/clusters/ ${NUM}.clusters
        cp cluster_sizes.png ${NUM}_cluster_sizes.png
        cp cluster_sizes.pdf ${NUM}_cluster_sizes.pdf
        mv ${NUM}_cluster_sizes.png ./workup3/logs/
        mv ${NUM}_cluster_sizes.pdf ./workup3/logs/
done
```

All clusters were combined using cat.

#### 6. Convert clusters into cool format

The conversion script will output the tab-delimited cluster file into cool binary format. The cluster file needs to be compressed first. The cluster resolution was size classes 2 to 100.
```
source activate cooler

gzip ./workup3/clusters/all_Dcor.clusters

/central/groups/Parker_lab/tools/sprite-pipeline/make_sprite_cooler.sh \
./workup3/clusters/all_Dcor.clusters chrom.sizes 2 100
```
Output is all_Dcor.clusters_2-100.norm_n2.1000.cool.

#### 7. Binning and Normalization in hicexplorer (https://hicexplorer.readthedocs.io/en/latest/)
Within hicexplorer, the cool formatted matrix was converted to h5 format. Different bin sizes were tested based on their tutorial for resolution of whole genome vs. single chromosome.
```
# get ready for hi-c map generation
source activate hicex

hicConvertFormat --matrices ./workup3/clusters/all_Dcor.clusters_2-100.norm_n2.1000.cool \
--outFileName ./workup3/clusters/all_Dcor.h5 --inputFormat cool --outputFormat h5 --load_raw_values

# merge different bin sizes
hicMergeMatrixBins -m ./workup3/clusters/all_Dcor.h5 -o ./workup3/clusters/all_Dcor_merged_nb50.h5 -nb 50
hicMergeMatrixBins -m ./workup3/clusters/all_Dcor.h5 -o ./workup3/clusters/all_Dcor_merged_nb10.h5 -nb 10
hicMergeMatrixBins -m ./workup3/clusters/all_Dcor.h5 -o ./workup3/clusters/all_Dcor_merged_nb80.h5 -nb 80
hicMergeMatrixBins -m ./workup3/clusters/all_Dcor.h5 -o ./workup3/clusters/all_Dcor_merged_nb30.h5 -nb 30
```

Create plot for removing low coverage and high coverage bins
```
hicCorrectMatrix diagnostic_plot -m ./workup3/clusters/all_Dcor_merged_nb50.h5 -o ./workup3/clusters/hic_diagnostic_nb50.png
#hicCorrectMatrix diagnostic_plot -m ./workup3/clusters/all_Dcor_merged_nb10.h5 -o ./workup3/clusters/hic_diagnostic_nb10.png
#hicCorrectMatrix diagnostic_plot -m ./workup3/clusters/all_Dcor_merged_nb80.h5 -o ./workup3/clusters/hic_diagnostic_nb80.png
#hicCorrectMatrix diagnostic_plot -m ./workup3/clusters/all_Dcor_merged_nb30.h5 -o ./workup3/clusters/hic_diagnostic_nb30.png
```

nb30:
![alt text](./hic_diagnostic_nb30.png "Diagnostic_plot")

Correct matrix with chosen thresholds. Run Knight-Ruiz matrix balancing algorithm (KR) or the iterative matrix correction (ICE).
```
hicCorrectMatrix correct -m ./workup3/clusters/all_Dcor_merged_nb30.h5 --filterThreshold -2 2 \
-o ./workup3/clusters/Dcor_hic_corrected_nb30.h5
```

#### 8. HiCAssembler (https://pypi.org/project/HiCAssembler/)
Replicated min scaffold length and bin sizes from their paper on two different Drosophila assemblies at different degrees of a fragmented starting assembly.

Before Assembly- after misassembly detection:

![alt text](./before_assembly_200k.png "Before_assembly")
```
/home/sak3097/.local/bin/assemble --matrix ./workup3/clusters/Dcor_hic_corrected_nb30.h5 \
--outFolder ./workup3/hi-c_assembly_nb30_scaff200k_it4 \
--fasta /central/groups/Parker_lab/final_genomes/Dalotia/Dcor_assembly_v1_190722.pruned_scaffolds.fa \
--min_scaffold_length 200000 --bin_size 10000 --num_processors 10 --misassembly_zscore_threshold -1.0 --num_iterations 4

/home/sak3097/.local/bin/assemble --matrix ./workup3/clusters/Dcor_hic_corrected_nb30.h5 \
--outFolder ./workup3/hi-c_assembly_nb30_scaff100k_it4 \
--fasta /central/groups/Parker_lab/final_genomes/Dalotia/Dcor_assembly_v1_190722.pruned_scaffolds.fa \
--min_scaffold_length 100000 --bin_size 5000 --num_processors 10 --misassembly_zscore_threshold -1.0 --num_iterations 4
```

I tried various combination of clusters (2-10, 2-50, 2-100, and 2-1000), matrix bin sizes (0, 10, 30, 50 and 80), iterations (2-5), minimum scaffold lengths (100000 to 300000) and assembly bin sizes (15000 to 5000). The clusters sizes had little impact and generally all contain large interchromosomal interactions compared to Hi-C based data. The matrix bin size and iterations were the tuning parameters to finalize the assembly, where moderately sized bins (30-50) and 4 iterations both gave the best final assemblies.

After Assembly (100k min scaffold size)

![alt text](./after_put_scaff_back.png "After_assembly")

Several of the final maps combined scaffolds 3-5 into one larger pseudomolecule. This was manually split based on the coordinates recovered from plotScaffoldInteractive from HiCAssembler.

Scaffold 4 was split into 3 and the remaining scaffolds were renamed based on their heatmap location.

```
bedtools getfasta -name -fi super_scaffolds.fa -bed split_scaff.bed -fo hic_scaff_split.fa
sed 's/\::.*$//' hic_scaff_split.fa > out.fa
mv out.fa hic_scaff_split.fa
```

##### 9. Gap-filling the SPRITE assembly with long-reads
There were 344 gaps and 465442 ambiguous bases.
```
# gapfilling with uncorrected nanopore reads
/central/groups/Parker_lab/tools/LR_Gapcloser_v1.1/LR_Gapcloser.sh -i hic_scaff_split.fa \
-l /central/groups/Parker_lab/raw_genome_data/1_Dalotia/minION/all_dalotia_minION.fasta \
-s n -t 10 -o gapFill_uncorHiC
```
After three iterations, there were 137 remaining gaps and 113,004 ambiguous bases.

### 10. Final round of polishing
Used Pilon:
```
#load correct python module
module load samtools/1.8

/central/groups/Parker_lab/tools/minimap2/minimap2 -t 10 -ax sr \
./gapFill_uncorHiC/iteration-3/gapclosed.fasta \
/central/scratchio/sak3097/beetle/1/filtered/1_filtered_R1.fastq \
/central/scratchio/sak3097/beetle/1/filtered/1_filtered_R2.fastq \
| samtools view -b - > /central/scratchio/sak3097/beetle/1/hic_nb30_100k_SR.bam

samtools sort -@ 10 -o /central/scratchio/sak3097/beetle/1/sorted_hic_nb30_100k_SR.bam \
/central/scratchio/sak3097/beetle/1/hic_nb30_100k_SR.bam

samtools index /central/scratchio/sak3097/beetle/1/sorted_hic_nb30_100k_SR.bam

java -Xmx240G -jar /central/groups/Parker_lab/tools/pilon-1.23.jar \
--genome ./gapFill_uncorHiC/iteration-3/gapclosed.fasta \
--frags /central/scratchio/sak3097/beetle/1/sorted_hic_nb30_100k_SR.bam \
--output Dcor_HiC_pilon1 --changes --diploid \
--threads 10 --chunksize 1000000 --tracks
```

169450181 reads mapped
164395802 properly paired
Mean total coverage: 173
17 gaps partially filled (134 gaps remain)
112,599 ambiguous bases

### 11. Re-run repeatModeler and repeatMasker
I repeated repeatModeler first to identify any possible new predictions based on the re-ordered and more contiguous genome assembly. The new predictions were combined with the previous repeat library and clustered with VSEARCH.
```
#combine repbase, mitetracker, and repeatmodeler
echo "## Combine repeat libraries"
cat /central/groups/Parker_lab/tools/RepeatMasker/util/repeatmasker.Coleoptera.fa \
/central/groups/Parker_lab/tools/MITE-Tracker/results/${NAME}_MT/families_nr.fasta.classified \
/central/groups/Parker_lab/SheilaK/sprite_dalotia/liftover/repeatMask_newModels/Dcor_SPRITE.repeatlib.fa \
/central/groups/Parker_lab/SheilaK/genome_assemblies/1_Dalotia/repeat_modeler/Dcor.repeatlib.fasta \
 > ${NAME}.ALL_2.fa

echo "$(grep -c '>' ${NAME}.ALL_2.fa) sequences in full library"

#cluster families
echo "## Cluster repeat libraries"
/central/groups/Parker_lab/tools/MITE-Tracker/vsearch-2.7.1/bin/vsearch \
--cluster_fast ${NAME}.ALL_2.fa --threads 4 --strand both \
--clusters ${NAME}_2.clust_temp --iddef 1 --id 0.8 --uc out.uc --centroids ${NAME}_2.repeatlib.fa

echo "$(grep -c '>' ${NAME}_2.repeatlib.fa) sequences in clustered library"
```
Repeatmasking increased from 11.58% to 13.51%. The most abundant TEs are LTR gypsy (1.21%), RC Helitron (2.3%) and DNA hAT (1.76%).

### 12. Liftover and gene annotation re-analysis
I used UCSC liftover utilties to transfer the previous annotation to the new genome assembly. The steps to create the liftover chain were taken from this bash script: https://hgwdev.gi.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt

627 genes were removed/split after the liftover process.

To re-identify these lost/split genes, the new assembly was run through GeMoMa analysis using the combined Nicophorous and Tribolium gene models and PASA using the same transcriptome as before.

Liftover, GeMoMa and PASA predictions were combined with evidenceModeler.

A final round of PASA was run to compare the evidenceModeler predicted genes against the transcriptome and add UTRs.

The 612 "new" gene models were blast against the respective database (NCBI, SwissProt, Cazy, MEROPS, EGGNOG) and combined with the previous annotation. Annotation was added to the GFF file using GAG.  
