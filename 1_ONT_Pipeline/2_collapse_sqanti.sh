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
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/iPSC_ABhinge
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_ONT_Pipeline/ipscABhinge_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/assist_ont_processing
export PATH=$PATH:${LOGEN_ROOT}/target_gene_annotation


##-------------------------------------------------------------------------

# merge alignment
echo "Collapsing..."
allfilteredmapped=($(ls ${WKD_ROOT}/3_align/*filtered.bam)) 
ls ${allfilteredmapped[@]}
source activate nanopore
samtools merge -f ${WKD_ROOT}/4_collapse/${samplename}_mapped.filtered.sorted.bam ${allfilteredmapped[@]}

# collapse
echo "Collapsing..."
echo "Output: ${WKD_ROOT}/4_collapse/${samplename}_collapsed.gff"
cd ${WKD_ROOT}/4_collapse
source activate isoseq3
isoseq3 collapse ${WKD_ROOT}/4_collapse/${samplename}_mapped.filtered.sorted.bam ${samplename}_collapsed.gff \
  --min-aln-coverage 0.85 --min-aln-identity 0.95 --do-not-collapse-extra-5exons \
  --log-level TRACE --log-file ${samplename}_collapsed.log

# demultiplex 
source activate nanopore
cd ${WKD_ROOT}/4_collapse
adapt_cupcake_to_ont.py ${WKD_ROOT}/3_align -o ${samplename} > adapt_cupcake.log

demux_cupcake_collapse.py \
  ${WKD_ROOT}/4_collapse/${samplename}_collapsed.read_stat.txt \
  ${WKD_ROOT}/3_align/${samplename}_sample_id.csv\
  --dataset=ont

# sqanti3
echo "Running SQANTI3..."
source activate sqanti2_py3
cd ${WKD_ROOT}/5_sqanti3
python $SQANTI3_DIR/sqanti3_qc.py ${WKD_ROOT}/4_collapse/${samplename}_collapsed.gff \
$GENOME_GTF $GENOME_FASTA -t 30 -fl ${WKD_ROOT}/4_collapse/demux_fl_count.csv \
--genename --skipORF --report skip &> ${samplename}_sqanti_qc.log

# sqanti3 filter 
filteringJson=$SQANTI3_DIR/utilities/filter/filter_default_reducecoverage.json
python $SQANTI3_DIR/sqanti3_filter.py rules ${samplename}_collapsed_classification.txt \
--faa=${samplename}_collapsed_corrected.fasta \
--gtf=${samplename}_collapsed_corrected.gtf \
-j=${filteringJson} --skip_report &> ${samplename}_sqanti_filter.log



# cpat
source activate nanopore
cd ${WKD_ROOT}/5_sqanti3
Rscript ${LOGEN_ROOT}/target_gene_annotation/identify_true_isoforms.R -s ${samplename}_collapsed_RulesFilter_result_classification.txt
seqtk subseq ipsc_collapsed_corrected.fasta ipsc_collapsed_RulesFilter_result_classification_isoform.txt > ipsc_collapsed.filtered.fasta
run_cpat ${WKD_ROOT}/5_sqanti3/${samplename}_collapsed.filtered.fasta ${samplename} ${WKD_ROOT}/6_cpat

# generate reference gtf for UNC13A and STMN2
grep -F -f target_genes.txt ipsc_collapsed_RulesFilter_result_classification_counts.txt  | awk '{print $1}'  > target_filtered_isoforms.txt
grep -F -f target_genes.txt $GENOME_GTF > gencode.v19.annotation_target.gtf
grep -F -f target_filtered_isoforms.txt ipsc_collapsed.filtered.gtf > target_ipsc_collapsed.filtered.gtf
grep -F -f ${WKD_ROOT}/5_sqanti3/target_filtered_isoforms.txt ${WKD_ROOT}/6_cpat/ipsc.ORF_prob.best.tsv > ${WKD_ROOT}/6_cpat/ipsc_target.ORF_prob.best.tsv
grep -F -f ${WKD_ROOT}/5_sqanti3/target_filtered_isoforms.txt ${WKD_ROOT}/6_cpat/ipsc.no_ORF.txt > ${WKD_ROOT}/6_cpat/ipsc_target.no_ORF.txt

# ficle
source activate sqanti2_py3
prefix=${WKD_ROOT}/5_sqanti3/ipsc_collapsed.filtered
gtfToGenePred ${prefix}.gtf ${prefix}.genePred
genePredToBed ${prefix}.genePred > ${prefix}.bed12
sort -k1,1 -k2,2n ${prefix}.bed12 > ${prefix}_sorted.bed12

ficle.py --gene=STMN2 \
    --reference=${WKD_ROOT}/7_ficle/STMN2_gencode.gtf \
    --input_bed=${WKD_ROOT}/5_sqanti3/${samplename}_collapsed.filtered_sorted.bed12 \
    --input_gtf=${WKD_ROOT}/5_sqanti3/${samplename}_collapsed.filtered.gtf  \
    ---input_class=${WKD_ROOT}/5_sqanti3/${samplename}_collapsed_RulesFilter_result_classification.txt \
    --cpat=<path/to/cpat_ORF_prob.best.tsv>  \
    --output_dir=<path/to/output/directory>