#!/bin/bash
#SBATCH --time=168:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=50G   # memory per CPU core
#SBATCH -J "cf_cyp"   # job name
#SBATCH --output "cf_cyp"

source /central/groups/Bi160/tools/anaconda/etc/profile.d/conda.sh
conda activate ete

dir=/central/groups/Parker_lab/SheilaK/alignments_biosynthesis
USER=sak3097
TEST=OG0001112
GENE=CYP

# branch-site tests
# gland meos
ete3 evol -t $dir/${TEST}_${GENE}/${TEST}_${GENE}_trimmed.contree --alg $dir/${TEST}_${GENE}/${TEST}_transcript_aligned.fas \
-o cyp_codeml_Nov22_results -v 3 --models bsA bsA1 --tests bsA,bsA1 \
--mark Ene_evm.model.scaffold18508.1,,Liom_evm.model.scaffold02141.1 --cpu 2
