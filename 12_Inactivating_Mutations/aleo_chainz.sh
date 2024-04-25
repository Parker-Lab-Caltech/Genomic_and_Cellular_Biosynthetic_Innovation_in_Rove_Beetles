#!/bin/bash
#SBATCH --time=48:00:00  # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=20G   # memory per CPU core
#SBATCH -J "ch31"   # job name
#SBATCH --output "ch31"

source /central/groups/Parker_lab/tools/anaconda/etc/profile.d/conda.sh
conda activate toga

#module load gcc/9.2.0
#module load R/4.0.3

Dcor_genome=/central/groups/Parker_lab/final_genomes/Dalotia/Dcor_assembly_v3_220711.combined.fasta

NUM=3
SPECIES=Ecitophya
query_genome=/central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/repeatMask/RM_output/${NUM}_scaffolds.reduced_1k.filled.masked.fa

./make_chains.py Dcor $SPECIES ${Dcor_genome} ${query_genome} --executor slurm --project_dir Dcor_${SPECIES}_Qhox --force_def \
--lastz /central/groups/Parker_lab/tools/lastz-distrib/bin/lastz --executor_queuesize 10 --executor_partition any \
--blastz_h 2000 --blastz_y 3400 --blastz_l 4000 --blastz_k 2200
