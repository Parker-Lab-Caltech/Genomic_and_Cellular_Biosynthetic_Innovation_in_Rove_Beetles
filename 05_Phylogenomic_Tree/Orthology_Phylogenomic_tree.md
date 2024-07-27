## **Orthologs and Phylogenomic Tree**

Sheila Kitchen, 5/28/2020

### Sequence preparation for OrthoFinder
1. Copy each species protein fasta file into a common directory (ie. ./prot_seqs).
2. For each species that will be added to the analysis, it is best to remove any trailing information in the fasta header after the gene id. Additionally, adding a taxa code to each species respective fasta file will help with downstream analyses. Modify the headers to remove extra information and add a taxa code (ex. >Geo_evm.model.scaffold00001.3)
3. Remove isoform variants in the genomes with cdhit:
```
/central/groups/Parker_lab/tools/cd-hit-v4.8.1-2019-0228/cd-hit -i ./prot_seqs/Agla.aa.fasta -o ./cd_hit_output/Agla_cdHIT_output.fa -c 0.98 -T 6 -aL 0.3 -aS 0.3
```

|Species| Number of Genes| Number Filtered\*|
|----------|-----:|---:|
|Agla|20,568|15,985|
|Apla|22,062|14,936|
|Atum|17,120|13,829|
|Dpon|19,975|13,715|
|Ldec|18,737|15,232|
|Nves|19,532|14,454|
|Otau|21,639|15,922|
|Tcas|22,549|15,636|

4. *Dalotia* predicted isoforms were removed based on the gene id. 15,818 remained out of 18,150 after filtering.

### Ortholog identification with OrthoFinder
Run OrthoFinder (v2.3.3):
```
conda activate orthofinder

orthofinder -t 16 -a 1 -M msa -S diamond -A mafft -T fasttree -f ./cd_hit_output
```
-M method for gene tree inference, msa= multiple sequence alignment
-A use MAFFT for the alignment
-T use fasttree for tree building
-t number of sequence search threads
-a number of orthofinder threads

With a new version of the *Dalotia* genome, we had to remove version 1 gene models with version 2 gene models.
To remove a species from the analysis point the analysis to the previous working directory where the species to remove has to be commented out of `SpeciesIDs.txt`. To add a new species or genome version indicate the directory location of the new taxa (in this case ./cd_hit_output/Dcor_v2):
```
orthofinder -b /central/groups/Parker_lab/SheilaK/orthoFind/cd_hit_output/OrthoFinder/Results_Oct13_1/WorkingDirectory -f ./cd_hit_output/Dcor_v2
```
Details of all the output from OrthoFinder can be found at https://github.com/davidemms/OrthoFinder

#### Single-copy Orthologs
Within the results of OrthoFinder, a list of single-copy orthologs is written to the `Orthogroups` directory. In our case, we identified only 276 single-copy orthologs. These were used to reconstruct the species tree.

#### Many-to-many orthologs
Additionally, the many-to-many orthologs can be parsed out of the Orthogroups.GeneCount.tsv within the `Orthogroups` directory. We set a threshold of at least 80% of the taxa to be present in an orthogroup to test the 1-1 (n =808) and 1-many (n=8012) orthologs for downstream analyses.

### Gather single-copy alignments into a new directory
```
dir=/central/groups/Parker_lab/SheilaK/orthoFind/

for GROUP in $(cat $dir/cd_hit_output/seq_list/Orthogroups_SingleCopyOrthologues.txt); do
        mkdir ./alignments/${GROUP}
        # copy alignment for OF into the new directory
        cp $dir/cd_hit_output/OrthoFinder/Results_Oct13_1/MultipleSequenceAlignments/${GROUP}.fa \
        $dir/alignments/${GROUP}
        # remove cdHIT_output
        sed 's/.*cdHIT_output_/>/' $dir/alignments/${GROUP}/${GROUP}.fa > $dir/alignments/${GROUP}/${GROUP}_rename.fa
        # remove Dcor_noIsoforms_aa_
        sed 's/.*Dcor_noIsoforms_aa_/>/' $dir/alignments/${GROUP}/${GROUP}_rename.fa > $dir/alignments/${GROUP}/${GROUP}_rename2.fa
        # extract fasta header to text file
        grep ">" $dir/alignments/${GROUP}/${GROUP}_rename2.fa > $dir/alignments/${GROUP}/${GROUP}_seq_head.txt
        # remove the caret symbol
        sed 's/>//' $dir/alignments/${GROUP}/${GROUP}_seq_head.txt > $dir/alignments/${GROUP}/${GROUP}_seq_head2.txt
        # pull out sequences from the transcriptomes/mRNA files
        fastafetch -f /central/groups/Parker_lab/SheilaK/orthoFind/nuc_seqs/all.nuc.fasta \
        -i /central/groups/Parker_lab/SheilaK/orthoFind/nuc_seqs/all.nuc.fasta.idx -F \
        -q $dir/alignments/${GROUP}/${GROUP}_seq_head2.txt \
        > $dir/alignments/${GROUP}/${GROUP}_transcripts.fasta
        # remove intermediate files
        rm $dir/alignments/${GROUP}/${GROUP}_seq_head.txt
        rm $dir/alignments/${GROUP}/${GROUP}_rename.fa
done
```

### Trim Single-copy Multiple Sequence Alignments
Trimal v 1.4.1 "tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment" (http://trimal.cgenomics.org/):

```
dir=/central/groups/Parker_lab/SheilaK/orthoFind/

for GROUP in $(cat $dir/cd_hit_output/seq_list/Orthogroups_SingleCopyOrthologues.txt); do
        #remove * (stop codon) from the alignment
        sed 's/[*]\+/-/g' $dir/alignments/${GROUP}/${GROUP}_rename2.fa \
        > $dir/alignments/${GROUP}/${GROUP}_rename3.fa
        #run trimal
        /central/groups/Parker_lab/tools/trimal-1.4.1/source/trimal \
        -in $dir/alignments/${GROUP}/${GROUP}_rename3.fa \
        -out $dir/alignments/${GROUP}/${GROUP}_trimmed.fas -keepheader -fasta -gappyout
        #remove everything but the taxa code from the alignment
        sed 's/\_.*$//' $dir/alignments/${GROUP}/${GROUP}_trimmed.fas > out.fas
        #rename the trimmed alignment
        mv out.fas $dir/alignments/${GROUP}/${GROUP}_trimmed.fas
        #copy trimmed fasta to a directory for the next step
        cp $dir/alignments/${GROUP}/${GROUP}_trimmed.fas \
        /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/fasconcat
done
```

### Combine the trimmed MSA into a single super-matrix
The multiple sequence alignments are combined using FasConCat v1.04 (https://github.com/PatrickKueck/FASconCAT-G)

```
cd /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/fasconcat

perl /central/groups/Parker_lab/tools/FASconCAT-G/FASconCAT-G_v1.04.pl -s -l

# convert the partition file to be compatible with MARE
sed 's/LG,/charset/' FcC_supermatrix_partition.txt > FcC_charset.txt
sed -i '/[^;] *$/s/$/;/' FcC_charset.txt
```

The FcC_supermatrix_partition.txt provides the order and positions of each orthogroup within the super-matrix. There 113802aa in 276 partitions/orthogroups.

### Matrix Reduction After Identifying Phylogenetic Informative Sites
MARE v0.1.2 (https://www.zfmk.de/en/research/research-centres-and-groups/mare)

```
#move into the working directory
cd /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/fasconcat

#run MARE
/central/groups/Parker_lab/tools/MARE_v0.1.2-rc/MARE FcC_charset.txt FcC_supermatrix.fas

# move into the output directory of MARE
cd ./results

#rename the output matrix to be read by FasConCat
mv FcC_supermatrix.fas_reduced FcC_supermatrix_reduced.fas

#convert formats to phylip for partitionfinder
perl /central/groups/Parker_lab/tools/FASconCAT-G/FASconCAT-G_v1.04.pl -s -o -p -p

# copy matrix and partition file into new directory to run partitionfinder
cp FcC_FcC_supermatrix_reduced.phy ../../partitionFinder/
cp FcC_charset.txt_reduced ../../partitionFinder/

```
After running MARE, the matrix was reduced to 178 orthogroups with 67,138 aa.

### Model selection by PartitionFinder v2.1.1
The best scheme (step 52) was identified using the settings below:
```
source activate funannotate

cd /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/partitionfinder

python /central/groups/Parker_lab/tools/partitionfinder-2.1.1/PartitionFinderProtein.py /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/partitionFinder --raxml --rcluster-max 1000 --rcluster-percent 10 -p 16
```

### Maximum Likelihood tree (IQtree)
The reduced matrix and best partition scheme are used to reconstruct the ML tree with 1000 bootstrap replicates.
```
cd /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/iqTree

# protein substitution model from PartitionFinder2 and ML tree construction
/central/groups/Parker_lab/tools/iqtree-1.6.8-Linux/bin/iqtree \
-s /central/groups/Parker_lab/SheilaK/orthoFind/single-copy_allTaxa/partitionFinder/FcC_FcC_supermatrix_reduced.phy \
-spp sc_partionFinder.nex -bb 1000 -nt 16 -pre sc_pf
```
The {sc_pf}.contree is the consensus tree of the 1000 bootstrap on the AIC and BIC scores.

### Coalescent-model ASTRAL on unrooted gene trees
"ASTRAL seeks to find the tree that maximizes the number of induced quartet trees in gene trees that are shared by the species tree" (https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md)

#### Construct unrooted gene trees for single-copy orthologs
This is constructed for the 276 orthogroups.
```
dir=/central/groups/Parker_lab/SheilaK/orthoFind/

for GROUP in $(cat $dir/cd_hit_output/seq_list/Orthogroups_SingleCopyOrthologues_sub.txt); do
       # ML tree construction
       cd $dir/alignments/${GROUP}
       /central/groups/Parker_lab/tools/iqtree-1.6.8-Linux/bin/iqtree -nt 6 -bb 1000 -s $dir/alignments/${GROUP}/${GROUP}_trimmed.fas -wbt -pre ${GROUP}
done
```
#### Combine all gene trees
```
dir=/central/groups/Parker_lab/SheilaK/orthoFind/

mkdir $dir/single-copy_allTaxa/astral
cd $dir/single-copy_allTaxa/astral

cat $dir/alignments/*/*.treefile > $dir/single-copy_allTaxa/astral/all_ML.tree
```

#### Run ASTRAL v5.6.3
```
java -jar /central/groups/Parker_lab/tools/Astral/astral.5.6.3.jar \
--input $dir/single-copy_allTaxa/astral/all_ML.tree \
--output $dir/single-copy_allTaxa/astral/all_ML_spp.tree 2> $dir/single-copy_allTaxa/astral/all_ML.log
```

Internal branch lengths are in coalescent units (discordance between the gene trees) and nodal support values are the local posterior probabilities.
