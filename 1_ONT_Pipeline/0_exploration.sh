#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/iPSC_ABhinge
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_ONT_Pipeline/ipscABinge_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing

# merge original mapped files (from Abhinge) to view on IGV
controlSam=(10833_E5_nes0_R1_pass.sorted.sam 10833_E5_nes0_R2_pass.sorted.sam 10833_E8_nes0_R1_pass.sorted.sam 10833_E8_nes0_R2_pass.sorted.sam)
alsSam=(10833_E5_nesM_R1_pass.sorted.sam 10833_E5_nesM_R2_pass.sorted.sam 10833_E8_nesM_R1_pass.sorted.sam 10833_E8_nesM_R2_pass.sorted.sam)

cd ${WKD_ROOT}/1_minimap/
#cat ${controlSam[@]} > control_pass.sorted.sam
#cat ${alsSam[@]} > als_pass.sorted.sam

source activate sqanti2_py3 
export PATH=$PATH:/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/software/cDNA_Cupcake/sequence/
sam_to_gff3.py ${WKD_ROOT}/1_minimap/control_pass.sorted.sam -s hg19
gffread control_pass.sorted.gff3 -T -o control_pass.sorted.gtf

sam_to_gff3.py ${WKD_ROOT}/1_minimap/als_pass.sorted.sam -s hg19
gffread als_pass.sorted.gff3 -T -o als_pass.sorted.gtf


#grep -w PB.29824 ${WKD_ROOT}/5_sqanti3/ipsc_collapsed_corrected.gtf > UNC13A_collapsed_corrected.gtf
#grep -w PB.16227 ${WKD_ROOT}/5_sqanti3/ipsc_collapsed_corrected.gtf > STMN2_collapsed_corrected.gtf
