#!/usr/bin/env Rscript


## ---------- functions -----------------

library("dplyr")

SC_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/iPSC_ABhinge"
LOGEN_ROOT = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/"
source(paste0(LOGEN_ROOT, "merge_characterise_dataset/run_ggtranscript.R"))
source(paste0(SC_ROOT,"/1_ONT_Pipeline/ipscABhinge_ont.config.R"))
cpat <- as.data.frame(fread("/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/6_cpat/ipsc_target.ORF_prob.best.tsv"))
colnames(cpat) <- c("ID","ORFname","tlength","ORF_strand","ORF1","ORF_length","ORF_frame","ORF_start","ORF_end ORF","Fickett Hexamer","Coding_prob")

## --- plot ggtranscript at a gene level ---

gene_ggtranscript <- function(gene,inputClassFile,inputTitle=NULL,inputCpat=NULL,cpatSpecies=NULL){
  
  Isodf <- data.frame(
    Isoform = unlist(Isodf <- list(
      Reference = unique(gtf$ref_target[gtf$ref_target$gene_name == gene & !is.na(gtf$ref_target$transcript_id), "transcript_id"]),
      ALS = as.character(as.character(inputClassFile[which(inputClassFile$associated_gene == gene & inputClassFile$Dataset == "als"),"isoform"])),
      Control = as.character(as.character(inputClassFile[which(inputClassFile$associated_gene == gene & inputClassFile$Dataset == "control"),"isoform"])),
      Both = as.character(as.character(inputClassFile[which(inputClassFile$associated_gene == gene & inputClassFile$Dataset == "Both"),"isoform"]))
    )),
    Category = rep(names(Isodf), lengths(Isodf))
  )
  Isodf$colour <- c(rep(NA,length(Isodf$Category[Isodf $Category != "DTE"])))
  #print(inputCpat)
  p <- ggTranPlots(inputgtf=gtf$merged, classfiles=inputClassFile,
                   isoList = c(as.character(Isodf$Isoform)),
                   selfDf = Isodf, gene=gene, inputCpat=inputCpat, cpatSpecies=cpatSpecies) + labs(title = inputTitle) 
  return(p)
}


gene_flCounts <- function(gene,inputClassFile,inputTitle=NULL){
  p <- inputClassFile %>% filter(associated_gene == gene) %>% 
    dplyr::select(isoform, contains("10833"),Dataset, structural_category) %>% 
    reshape2::melt(variable.name = "sample", value.name = "reads", id = c("isoform","Dataset","structural_category")) %>% 
    mutate(id = paste0("LR.",gene,".",word(isoform,c(3),sep=fixed(".")))) %>%
    ggplot(., aes(x = reorder(id, -reads), y = reads)) + geom_boxplot(aes(colour = structural_category)) + 
    facet_grid(~Dataset, scales = "free", space = "free") + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    labs(x = "Transcript", y = "Number of full-length ONT reads",title=inputTitle) +
    theme( strip.background = element_blank(), legend.position = "top")
  
  return(p)
}


## ---------- analysis -----------------

# Tracks for each target gene
pGeneTranscripts <- list()
pGeneCounts <- list()
for(gene in c("UNC13A","CELF5","ELAVL3","CACNA1E","KCNQ2")){
  pGeneTranscripts[[gene]] <- gene_ggtranscript(gene,class.files,inputTitle=gene) 
  pGeneCounts[[gene]] <- gene_flCounts(gene,class.files,inputTitle=gene) 
}
pGeneTranscripts$UNC13A <- pGeneTranscripts$UNC13A + annotate("rect", xmin = 12000, xmax = 12200, ymin = -1, ymax = 11, alpha = 0.2, fill = "green")


# --- STMN2 ---
# plot of coding potential by transcript length for STMN2 isoforms 
# to observe if there is a correlation between the short isoforms and coding potential
cpat[cpat$ID %in% class.files[class.files$associated_gene == "STMN2","isoform"],] %>% 
  left_join(., class.files[,c("isoform","length","structural_category")], by = c("ID" = "isoform")) %>% 
  ggplot(., aes(x = Coding_prob, y = length)) + geom_point(aes(colour = structural_category)) +
  theme_classic() + labs(x = "Coding probability", y = "Transcript length (bp)") +
  geom_vline(xintercept = 0.364, linetype="dotted")

# filter for non-coding transcripts according to cpat (default threshold for human = 0.364)
STMN2np =  c(cpat[cpat$ID %in% class.files[class.files$associated_gene == "STMN2","isoform"] & cpat$Coding_prob < 0.364,"ID"],"PB.16234.4701")
class.files.STMN2np = class.files[class.files$isoform %in% STMN2np, ]
pGeneTranscripts$STMN2 <- gene_ggtranscript("STMN2",class.files.STMN2np,"STMN2") 
pGeneCounts$STMN2 <- gene_flCounts("STMN2",class.files.STMN2np,inputTitle="STMN2") 
pGeneTranscripts$STMN2CE <- gene_ggtranscript("STMN2",class.files[class.files$isoform %in% SMN2CEIso, ],"STMN2",inputCpat=cpat,cpatSpecies="human") 

nrow(class.files[class.files$associated_gene == "STMN2",])
nrow(class.files[class.files$associated_gene == "STMN2" & class.files$structural_category == "FSM",])
# number of isoforms with crytpic exon with start coordinate "80529057"
SMN2CE <- read.table("/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/5_sqanti3/SMN2_cryptic_exon.gtf")
SMN2NE <- read.table("/lustre/projects/Research_Project-MRC148213/sl693/ipscABhinge/5_sqanti3/SMN2_1stnovel_exon.gtf")
SMN2CEIso <- unique(SMN2CE$V13)
SMN2NEIso <- unique(SMN2NE$V13)
length(SMN2CEIso) 
length(SMN2NEIso)
unique(SMN2NE$V4 - SMN2NE$V5)
setdiff(SMN2NEIso,SMN2CEIso)
class.files[class.files$isoform %in% unique(SMN2CE$V13),c("isoform","control_sum_FL","als_sum_FL")]
class.files[class.files$isoform %in% unique(SMN2NE$V13),c("isoform","control_sum_FL","als_sum_FL")]

setdiff(SMN2CEIso,STMN2np)
gene_flCounts("STMN2",class.files[class.files$isoform %in% c("PB.16234.4908"),])

STMN2Fl <- class.files[class.files$isoform %in% c("PB.16234.65","PB.16234.132","PB.16234.219","PB.16234.973","PB.16234.1823"),] %>%
  dplyr::select(isoform, contains("10833"),structural_category) %>% 
  reshape2::melt(variable.name = "sample", value.name = "reads", id = c("isoform","structural_category")) %>% 
  mutate(Dataset = word(sample,c(-1),sep=fixed("_"))) %>% 
  mutate(id = paste0("LR.",gene,".",word(isoform,c(3),sep=fixed(".")))) %>%
  ggplot(., aes(x = reorder(id, -reads), y = reads)) + geom_boxplot(aes(colour = structural_category)) + geom_point() +
  facet_grid(~Dataset, scales = "free", space = "free") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Transcript", y = "Number of full-length ONT reads") +
  theme( strip.background = element_blank(), legend.position = "top")

SMN2AF <- paste0("PB.16234.",c("4847","4850","4851","4852","4853","4854","4855","4857","4858","4861","4863"))
class.files[class.files$isoform %in% SMN2AF,c("isoform","control_sum_FL","als_sum_FL")]

## ---------- output -----------------

pdf(paste0(SC_ROOT, "/Figs/Tracks_UNC13A_CELF5_ELAVL3_KNCQ2.pdf"), width=10, height=8)
for(n in c(1,2,3,5)){print(pGeneTranscripts[n])}
dev.off()

pdf(paste0(SC_ROOT, "/Figs/Tracks_CACNA1E_STMN2.pdf"), width=15, height=15)
for(n in c(4,6)){print(pGeneTranscripts[n])}
dev.off()

pdf(paste0(SC_ROOT, "/Figs/transcriptFLCounts.pdf"),width=10)
pGeneCounts 
dev.off()

png(paste0(SC_ROOT, "/Figs/STMN2FLCounts.png"))
STMN2Fl 
dev.off()

png(paste0(SC_ROOT, "/Figs/Tracks_STMN2.png"), width = 600, height = 900)
pGeneTranscripts$STMN2CE
dev.off()
