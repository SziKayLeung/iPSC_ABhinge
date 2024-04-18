#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working WKD_ROOTectory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=5:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --mem=200G # specify bytes memory to reserve



# Batch2: realign with pbmm2 align and filter alignment 

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/iPSC_ABhinge
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_ONT_Pipeline/ipscABhinge_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

##-------------------------------------------------------------------------

source activate nanopore
cd ${WKD_ROOT}/5_sqanti3
Rscript ${LOGEN_ROOT}/target_gene_annotation/identify_true_isoforms.R -s ${samplename}_collapsed_RulesFilter_result_classification.txt
seqtk subseq ipsc_collapsed_corrected.fasta ipsc_collapsed_RulesFilter_result_classification_isoform.txt > ipsc_collapsed.filtered.fasta
run_cpat ${WKD_ROOT}/5_sqanti3/${samplename}_collapsed.filtered.fasta ${samplename} ${WKD_ROOT}/6_cpat
