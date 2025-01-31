#!/usr/bin/env Rscript
## ----------Script-----------------
##
## Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
## List of genes that are detected solely in the control vs als samples
## qPCR targets for Barry
## --------------------------------

## ---------- packages -----------------

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))

## ------------------------------------

dir <- "/lustre/projects/Research_Project-MRC148213/lsl693/ipscABhinge/5_sqanti3"

# read classification file
class.file <- fread(paste0(dir, "ipsc_collapsed_RulesFilter_result_classification_counts.txt"), data.table = F)
class.file <- class.file %>% filter(filter_result == "Isoform")

# tabulate the number of reads per gene 
geneCounts <- class.file %>% select( associated_gene, contains("mapped")) %>%
  group_by(associated_gene) %>%
  summarise_each(list(sum)) %>% 
  tibble::column_to_rownames(., var = "associated_gene")

# for each gene, tabulate the number of total reads in control vs als samples
controlGeneCounts <- as.data.frame(geneCounts %>% select(contains("control")) %>% apply(.,1,sum))
alsGeneCounts <- as.data.frame(geneCounts %>% select(contains("als")) %>% apply(.,1,sum))
geneCountGroup <- cbind(controlGeneCounts, alsGeneCounts)
colnames(geneCountGroup) <- c("control","als")

# create a dataset column 
identify_dataset_by_counts <- function(col1,col2,name1,name2){
  
  col1 = as.numeric(col1)
  col2 = as.numeric(col2)
  
  if(col1 > 0 & col2 > 0){return("Both")
  }else if(col1 == 0 & col2 > 0){return(name2)
  }else if(col1 > 0 & col2 == 0){return(name1)
  }else{return("NA")}
  
}

geneCountGroup$dataset <- apply(geneCountGroup, 1, function(x) 
  identify_dataset_by_counts (x[["control"]], x[["als"]], "Control","ALS"))

## ------------------------------------

# write output
geneCountGroup <- geneCountGroup %>% tibble::rownames_to_column(., var = "associated_gene")
write.table(geneCountGroup,paste0(dir, "ipsc_geneFLReadCounts.txt"), quote = F, sep = "\t", row.names = F)