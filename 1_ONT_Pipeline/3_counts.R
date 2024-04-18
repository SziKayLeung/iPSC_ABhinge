
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "/compare_datasets/dataset_identifer.R"))

counts = data.table::fread("/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/4_collapse/demux_fl_count.csv")
phenotype <- read.table(paste0(rootDir,"0_metadata/sample_metadata.txt"),sep="\t",header=T)
phenotype <- phenotype %>% mutate(col = paste0(sampleID,"_",phenotype))


names(counts)[-1] <- as.character(phenotype$col[match(names(counts), phenotype$sampleID)])[-1]
names(counts)[1] <- "isoform"

for(i in unique(phenotype$phenotype)){
  message("Summing reads for phenotype: ", i)
  colname = paste0(i,"_sum_FL")
  counts[[colname]] <- apply(counts %>% select(contains(i)),1,sum)
}

dataset1=as.character(unique(phenotype$phenotype)[1])
dataset2=as.character(unique(phenotype$phenotype)[2])

rootDir <- "/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/"
class.names.files <- paste0(rootDir,"5_sqanti3/ipsc_collapsed_RulesFilter_result_classification.txt")
class.files <- as.data.frame(data.table::fread(class.names.files))
class.files <- class.files[class.files$filter_result == "Isoform",]

counts$dataset <- apply(counts, 1, function(x) identify_dataset_by_counts(x[[paste0(dataset1,"_sum_FL")]], x[[paste0(dataset2,"_sum_FL")]], dataset1,dataset2))
class.files <- class.files %>% dplyr::select(-contains("FL."))
class.files <- merge(class.files, counts, by = "isoform", all.x = T)
write.table(class.files, "/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/5_sqanti3/ipsc_collapsed_RulesFilter_result_classification_counts.txt",sep="\t",quote=F,row.names=F)
