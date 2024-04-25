#!/bin/bash
#SBATCH --time 120:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8G   # memory per CPU core
#SBATCH -J "bayes_mcmcR2"   # job name
#SBATCH --output "bayes_mcmcR2"

source activate funannotate

#mcmctree, calculate Gradient and Hessian to Approximate the Likelihood
#/central/groups/Parker_lab/tools/paml4.9i/src/mcmctree mcmctree_outBV.ctl

#/central/groups/Parker_lab/tools/paml4.9i/src/codeml tmp0001.ctl

# run 2
/central/groups/Parker_lab/tools/paml4.9i/src/mcmctree mcmctree.ctl
