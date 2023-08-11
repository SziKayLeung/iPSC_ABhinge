#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working WKD_ROOTectory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mem=200G # specify bytes memory to reserve
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_collapse_sqanti3.o
#SBATCH --error=2_collapse_sqanti3.e


# Batch2: realign with pbmm2 align and filter alignment 

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/RWest
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_ONT_Pipeline/rwest_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing


##-------------------------------------------------------------------------

# merge alignment
echo "Collapsing..."
allfilteredmapped=($(ls ${WKD_ROOT}/5_align/*filtered.bam)) 
ls ${allfilteredmapped[@]}
source activate nanopore
samtools merge -f ${WKD_ROOT}/6_collapse/${samplename}_mapped.filtered.sorted.bam ${allfilteredmapped[@]}

# collapse
echo "Collapsing..."
echo "Output: ${WKD_ROOT}/6_collapse/${samplename}_collapsed.gff"
cd ${WKD_ROOT}/6_collapse
source activate isoseq3
isoseq3 collapse ${WKD_ROOT}/6_collapse/${samplename}_mapped.filtered.sorted.bam ${samplename}_collapsed.gff \
  --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
  --log-level TRACE --log-file ${samplename}_collapsed.log

# demultiplex 
source activate nanopore
cd ${WKD_ROOT}/6_collapse
adapt_cupcake_to_ont.py ${WKD_ROOT}/5_align -o ${samplename} > adapt_cupcake.log

demux_cupcake_collapse.py \
  ${WKD_ROOT}/6_collapse/${samplename}_collapsed.read_stat.txt \
  ${WKD_ROOT}/5_align/${samplename}_sample_id.csv\
  --dataset=ont

# sqanti3
echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${WKD_ROOT}/7_sqanti3
python $SQANTI3_DIR/sqanti3_qc.py ${WKD_ROOT}/6_collapse/${samplename}_collapsed.gff \
$GENOME_GTF $GENOME_FASTA -t 30 -fl ${WKD_ROOT}/6_collapse/demux_fl_count.csv \
--genename --skipORF --report skip &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
--faa=${samplename}_collapsed_corrected.fasta \
--gtf=${samplename}_collapsed_corrected.gtf \
-j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log