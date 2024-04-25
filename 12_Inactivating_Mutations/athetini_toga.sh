#!/bin/bash
#SBATCH --time=48:00:00  # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=20G   # memory per CPU core
#SBATCH -J "at_TOGA41"   # job name
#SBATCH --output "at_TOGA41"

source /central/groups/Parker_lab/tools/anaconda/etc/profile.d/conda.sh
conda activate toga

module load gcc/9.2.0
#module load R/4.0.3

Dcor_genome=/central/groups/Parker_lab/final_genomes/Dalotia/Dcor_assembly_v3_220711.combined.fasta

NUM=3
SPECIES=Ecitophya
query_genome=/central/groups/Parker_lab/SheilaK/genome_assemblies/${NUM}_${SPECIES}/repeatMask/RM_output/${NUM}_scaffolds.reduced_1k.filled.masked.fa


#./toga.py ../make_lastz_chains/Dcor_${SPECIES}_chain/Dcor.${SPECIES}.allfilled.chain.gz \
# /central/groups/Parker_lab/final_genomes/Dalotia/Dcor_agat.bed \
#/central/groups/Parker_lab/tools/make_lastz_chains/Dcor_${SPECIES}_chain/Dcor.2bit \
#/central/groups/Parker_lab/tools/make_lastz_chains/Dcor_${SPECIES}_chain/${SPECIES}.2bit \
#--kt --pn dcor_${SPECIES}_TOGA_frag --nc /central/groups/Parker_lab/tools/TOGA/nextflow_config_files \
#--cb 8,16,32,64,128 --chn 10 --cjn 300 --ms --fragmented_genome

./toga.py ../make_lastz_chains/Dcor_${SPECIES}_Qhox/Dcor.${SPECIES}.allfilled.chain.gz \
 /central/groups/Parker_lab/final_genomes/Dalotia/Dcor_agat.bed \
/central/groups/Parker_lab/tools/make_lastz_chains/Dcor_${SPECIES}_Qhox/Dcor.2bit \
/central/groups/Parker_lab/tools/make_lastz_chains/Dcor_${SPECIES}_Qhox/${SPECIES}.2bit \
--kt --pn dcor_${SPECIES}_TOGA_Qhox --nc /central/groups/Parker_lab/tools/TOGA/nextflow_config_files \
--cb 8,16,32,64,128 --chn 10 --cjn 300 --ms --fragmented_genome
