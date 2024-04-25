#!/bin/bash
#SBATCH --time=120:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "sel_cyp_fg"   # job name
#SBATCH --output "sel_cyp_fg"

source /central/groups/Bi160/tools/anaconda/etc/profile.d/conda.sh
conda activate bi160_hyphy

dir=/central/groups/Parker_lab/SheilaK/alignments_biosynthesis
USER=sak3097
TEST=OG0001112
GENE=CYP

# selected branch test
HYPHYMPI absrel \
--alignment $dir/${TEST}_${GENE}/${TEST}_transcript_aligned.fas \
--tree $dir/${TEST}_${GENE}/selection_tests/absrel_selected_branches/cyp4g_selected.nwk --branches Foreground \
--output $dir/${TEST}_${GENE}/selection_tests/absrel_selected_branches/selected_cyp4g.ABSREL.json CPU=8
