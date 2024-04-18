#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working WKD_ROOTectory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=3_rerun_sqanti.o
#SBATCH --error=3_rerun_sqanti.e


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

# sqanti3
echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${WKD_ROOT}/5_sqanti3
# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
--faa=${samplename}_collapsed_corrected.fasta \
--gtf=${samplename}_collapsed_corrected.gtf \
-j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log

source activate nanopore
cd ${WKD_ROOT}/5_sqanti3
Rscript ${LOGEN_ROOT}/target_gene_annotation/identify_true_isoforms.R -s ${samplename}_collapsed_RulesFilter_result_classification.txt
seqtk subseq ipsc_collapsed_corrected.fasta ipsc_collapsed_RulesFilter_result_classification_isoform.txt > ipsc_collapsed.filtered.fasta
run_cpat ${WKD_ROOT}/5_sqanti3/${samplename}_collapsed.filtered.fasta ${samplename} ${WKD_ROOT}/6_cpat


# ficle
source activate sqanti2_py3
prefix=${WKD_ROOT}/5_sqanti3/ipsc_collapsed.filtered
gtfToGenePred ${prefix}.gtf ${prefix}.genePred
genePredToBed ${prefix}.genePred > ${prefix}.bed12
sort -k1,1 -k2,2n ${prefix}.bed12 > ${prefix}_sorted.bed12

# SMN2
grep -w 80529057 ipsc_collapsed.filtered.gtf > SMN2_cryptic_exon.gtf
grep -w 80526058 ipsc_collapsed.filtered.gtf > SMN2_1stnovel_exon.gtf
