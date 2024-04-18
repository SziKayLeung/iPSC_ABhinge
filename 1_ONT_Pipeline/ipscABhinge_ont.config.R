
## ---------- packages -----------------

LOGEN <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen"
source(paste0(LOGEN,"/transcriptome_stats/read_sq_classification.R"))
source(paste0(LOGEN,"/target_gene_annotation/summarise_gene_stats.R"))
source(paste0(LOGEN,"/compare_datasets/dataset_identifer.R"))


## ---------- input -----------------

# directory names
rootDir <- "/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/"


# input phenotype
phenotype <- read.table(paste0(rootDir,"0_metadata/sample_metadata.txt"),sep="\t")

# SQANTI classification files
#unfiltered.class.names.files <- paste0(rootDir,"5_sqanti3/ipsc_collapsed_classification.txt")
#unfiltered.class.files <- read.table(unfiltered.class.names.files, sep = "\t", as.is = T, header = T)

# list of isoforms
#isoformsFiltered <- read.table(paste0(rootDir,"5_sqanti3/ipsc_collapsed_RulesFilter_result_classification_isoform.txt"))
#class.files <- unfiltered.class.files[unfiltered.class.files$isoform %in% isoformsFiltered$V1,]

class.names.files <- paste0(rootDir,"5_sqanti3/ipsc_collapsed_RulesFilter_result_classification_counts.txt")
class.files <- SQANTI_class_preparation(class.names.files,standard="all")


# gtf
gtf <- list(
  filtered = rtracklayer::import(paste0(rootDir,"5_sqanti3/target_ipsc_collapsed.filtered.gtf")),
  ref_target = rtracklayer::import(paste0(rootDir,"5_sqanti3/gencode.v19.annotation_target.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))
gtf$merged <- rbind(gtf$filtered[,c("seqnames","strand","start","end","type","transcript_id","gene_id")] ,
                    gtf$ref_target[,c("seqnames","strand","start","end","type","transcript_id","gene_id")])
