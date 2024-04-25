#!/bin/bash
#SBATCH --time=120:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "selCYP"   # job name
#SBATCH --output "selCYP"

# load modules and activate conda environment
module load samtools/1.8
source /central/groups/Bi160/tools/anaconda/etc/profile.d/conda.sh
conda activate bi160_hyphy

# set variables
dir=/central/groups/Parker_lab/SheilaK/alignments_biosynthesis
USER=sak3097
TEST=OG0001112
GENE=CYP

# align protein sequences
mafft --thread 6 --amino \
--auto $dir/${TEST}_${GENE}/sub_seqs/${TEST}_unaligned.aa.fa \
 > $dir/${TEST}_${GENE}/sub_seqs/${TEST}_aligned.aa.fas

#reorder mRNA sequences
grep ">" $dir/${TEST}_${GENE}/sub_seqs/${TEST}_aligned.aa.fas > $dir/${TEST}_${GENE}/sub_seqs/$gene.list
sed 's/>//'  $dir/${TEST}_${GENE}/sub_seqs/$gene.list >  $dir/${TEST}_${GENE}/sub_seqs/$gene.list2

# change multiline to single line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $dir/${TEST}_${GENE}/sub_seqs/${TEST}_unaligned.CDS.fa > \
$dir/${TEST}_${GENE}/sub_seqs/${TEST}_unaligned_single.CDS.fa

# index file
samtools faidx $dir/${TEST}_${GENE}/sub_seqs/${TEST}_unaligned_single.CDS.fa

samtools faidx $dir/${TEST}_${GENE}/sub_seqs/${TEST}_unaligned_single.CDS.fa $(cat  $dir/${TEST}_${GENE}/sub_seqs/$gene.list2) \
> $dir/${TEST}_${GENE}/sub_seqs/${TEST}_reordered.cds.fasta

# codon align mRNA sequences using the protein alignment
tranalign $dir/${TEST}_${GENE}/sub_seqs/${TEST}_reordered.cds.fasta \
$dir/${TEST}_${GENE}/sub_seqs/${TEST}_aligned.aa.fas \
-outseq $dir/${TEST}_${GENE}/sub_seqs/${TEST}_transcript_aligned.fas

# trim protein alignment
#trimal -in $dir/${TEST}_${GENE}/sub_seqs/${TEST}_aligned.aa.fas \
#-out $dir/${TEST}_${GENE}/sub_seqs/${TEST}_trimmed.fas -keepheader -fasta -gappyout

# calculate ML trees for protein alignment
#cd $dir/${TEST}_${GENE}/sub_seqs/

#iqtree -nt 6 -bb 1000 -s $dir/${TEST}_${GENE}/${TEST}_aligned.aa.fas -wbt -mset LG,WAG,JTT,Dayhoff,Q.insect -pre ${TEST}_${GENE}
#iqtree -nt 6 -bb 1000 -s $dir/${TEST}_${GENE}/sub_seqs/${TEST}_trimmed.fas -wbt -mset LG,WAG,JTT,Dayhoff,Q.insect -pre ${TEST}_${GENE}_trimmed_042123
