#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=4:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-1 # 2 barcodes
#SBATCH --output=1_cutadapt_minimap2_tclean-%A_%a.o
#SBATCH --error=1_cutadapt_minimap2_tclean-%A_%a.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/iPSC_ABhinge
source $SC_ROOT/1_ONT_Pipeline/ipscABinge_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh

sample=${SAMPLE_NAMES[${SLURM_ARRAY_TASK_ID}]}

##-------------------------------------------------------------------------

# convert sam file to bam file
source activate nanopore
samtools view -h ${WKD_ROOT}/1_minimap/${sample}_pass.sorted.bam > ${WKD_ROOT}/1_minimap/${sample}_pass.sorted.sam

# run transcript clean on aligned reads
run_transcriptclean ${WKD_ROOT}/1_minimap/${sample}_pass.sorted.sam ${WKD_ROOT}/2_tclean

# re-align with pbmm2
run_pbmm2 ${WKD_ROOT}/2_tclean/${sample}/${sample}_clean.fa ${WKD_ROOT}/3_align

# filter_alignment <input_name> <input_mapped_dir>
# output = ${sample}_mapped.filtered.bam, ${sample}_mapped.filtered.sorted.bam
filter_alignment ${sample}_mapped ${WKD_ROOT}/3_align