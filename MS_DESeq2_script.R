#DESeq2 analysis for microbiota data.


locols <- c('#FDF8FF' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
hicols <- c("#FFFCF3" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")



##### Also running ANCOM-BC in qiime2. #####

##### Make sure you have the appropriate unrarefied physeq objects #####
#Starting point is psexp, which is all experimental samples unrarefied

#   ##### Example subsetting for differential abundnace analysis #####
pspre1lotreat.nr <- subset_samples(psexp, day %in% c("D-5", "D1") &  microb=="LO" & treatment=="worm")
pspre1lotreat.nr <- prune_taxa(taxa_sums(pspre1lotreat.nr) >= 1, pspre1lotreat.nr)

pspre1loctrl.nr <- subset_samples(psexp, day %in% c("D-5", "D1") &  microb=="LO" & treatment=="control")
pspre1loctrl.nr <- prune_taxa(taxa_sums(pspre1loctrl.nr) >= 1, pspre1loctrl.nr)

pspre1hitreat.nr <- subset_samples(psexp, day %in% c("D-5", "D1") &  microb=="HI" & treatment=="worm")
pspre1hitreat.nr <- prune_taxa(taxa_sums(pspre1hitreat.nr) >= 1, pspre1hitreat.nr)

pspre1hictrl.nr <- subset_samples(psexp, day %in% c("D-5", "D1") &  microb=="HI" & treatment=="control")
pspre1hictrl.nr <- prune_taxa(taxa_sums(pspre1hictrl.nr) >= 1, pspre1hictrl.nr)
#

##### DESeq2 #####
#BiocManager::install(version='devel') ##-Updates BiocManager if need be.

library(phyloseq)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("DESeq2")
packageVersion("DESeq2") #1.34.0 for pre41, and 1.39.8 for pre14, 1.44.0 for all adjacent timepoints. 
#BiocManager::install("microbiome") 

library(DESeq2)
library(microbiome)
library(ashr)
library(apeglm)

tax_table <- phyloseq::tax_table

#Tutorial on DESeq2 from phyloseq.
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

#If you want to aggregate the taxa at a certain level, you can easily do so.
phy_fam <- tax_glom(pspre14loctrl.nr, taxrank="Family", NArm=TRUE) #outdated command was aggregate_taxa
#
###### Start DESeq #####
#   ##### pretreatment to D1 #####

#Run DESeq on lotreat
dsps = phyloseq_to_deseq2(pspre1lotreat.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D1_vs_D.5", type="ashr") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre1lotreat.nr)[rownames(resOrdered), ], "matrix"))
sigtablotreat = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1 ), ]
sigtablotreat = cbind(as(sigtablotreat, "data.frame"), as(tax_table(pspre1lotreat.nr)[rownames(sigtablotreat), ], "matrix"))
sigtablotreat$name <- paste(sigtablotreat$Family, sep=" ", rownames(sigtablotreat))
dim(sigtablotreat)

#Run DESeq on hitreat
dsps = phyloseq_to_deseq2(pspre1hitreat.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D1_vs_D.5", type="ashr") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre1hitreat.nr)[rownames(resOrdered), ], "matrix"))
sigtabhitreat = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1 ), ]
sigtabhitreat = cbind(as(sigtabhitreat, "data.frame"), as(tax_table(pspre1hitreat.nr)[rownames(sigtabhitreat), ], "matrix"))
sigtabhitreat$name <- paste(sigtabhitreat$Family, sep=" ", rownames(sigtabhitreat))
dim(sigtabhitreat)

#Run DESeq on loctrl
dsps = phyloseq_to_deseq2(pspre1loctrl.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D1_vs_D.5", type="ashr") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre1loctrl.nr)[rownames(resOrdered), ], "matrix"))
sigtabloctrl = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1 ), ]
sigtabloctrl = cbind(as(sigtabloctrl, "data.frame"), as(tax_table(pspre1loctrl.nr)[rownames(sigtabloctrl), ], "matrix"))
sigtabloctrl$name <- paste(sigtabloctrl$Family, sep=" ", rownames(sigtabloctrl))
dim(sigtabloctrl)

#Run DESeq on hictrl
dsps = phyloseq_to_deseq2(pspre1hictrl.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D1_vs_D.5", type="ashr") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre1hictrl.nr)[rownames(resOrdered), ], "matrix"))
sigtabhictrl = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1 ), ]
sigtabhictrl = cbind(as(sigtabhictrl, "data.frame"), as(tax_table(pspre1hictrl.nr)[rownames(sigtabhictrl), ], "matrix"))
sigtabhictrl$name <- paste(sigtabhictrl$Family, sep=" ", rownames(sigtabhictrl))
dim(sigtabhictrl)

#      ##### Bind all pretreatment to D1 results together #####

#No significant results. 

#   ##### pretreatment to D10 #####

#Run DESeq on lotreat
dsps = phyloseq_to_deseq2(pspre10lotreat.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D10_vs_D.5", type="ashr") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre10lotreat.nr)[rownames(resOrdered), ], "matrix"))
sigtablotreat = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1 ), ]
sigtablotreat = cbind(as(sigtablotreat, "data.frame"), as(tax_table(pspre10lotreat.nr)[rownames(sigtablotreat), ], "matrix"))
sigtablotreat$name <- paste(sigtablotreat$Family, sep=" ", rownames(sigtablotreat))
sigtablotreat$MT <- "LO_worm"
sigtablotreat$comparison <- "pre10"
dim(sigtablotreat)

#Run DESeq on hitreat
dsps = phyloseq_to_deseq2(pspre10hitreat.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D10_vs_D.5", type="ashr") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre10hitreat.nr)[rownames(resOrdered), ], "matrix"))
sigtabhitreat = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1  ), ]
sigtabhitreat = cbind(as(sigtabhitreat, "data.frame"), as(tax_table(pspre10hitreat.nr)[rownames(sigtabhitreat), ], "matrix"))
sigtabhitreat$name <- paste(sigtabhitreat$Family, sep=" ", rownames(sigtabhitreat))
sigtabhitreat$MT <- "HI_worm"
sigtabhitreat$comparison <- "pre10"
dim(sigtabhitreat)

#Run DESeq on loctrl
dsps = phyloseq_to_deseq2(pspre10loctrl.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D10_vs_D.5", type="apeglm") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(pspre10loctrl.nr)[rownames(resOrdered), ], "matrix"))
sigtabloctrl = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1), ]
sigtabloctrl = cbind(as(sigtabloctrl, "data.frame"), as(tax_table(pspre10loctrl.nr)[rownames(sigtabloctrl), ], "matrix"))
sigtabloctrl$name <- paste(sigtabloctrl$Family, sep=" ", rownames(sigtabloctrl))
sigtabloctrl$MT <- "LO_sham"
sigtabloctrl$comparison <- "pre10"
dim(sigtabloctrl)

#Run DESeq on hictrl
dsps = phyloseq_to_deseq2(pspre10hictrl.nr, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D10_vs_D.5", type="apeglm") #change coef
summary(resLFC)
resOrdered <- resLFC[order(resLFC$log2FoldChange),]
sigtabhictrl = resOrdered[which(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) >= 1 ), ]
sigtabhictrl = cbind(as(sigtabhictrl, "data.frame"), as(tax_table(pspre10hictrl.nr)[rownames(sigtabhictrl), ], "matrix"))
sigtabhictrl$name <- paste(sigtabhictrl$Family, sep=" ", rownames(sigtabhictrl))
sigtabhictrl$MT <- "HI_sham"
sigtabhictrl$comparison <- "pre10"
dim(sigtabhictrl)

##### Merge outputs from ANCOM-BC and DESeq2 #####

library(tidyverse)
library(stringr)
library(writexl)

#Outside of R, I renamed all files so that they had their parent file name attached. 
#Then I pulled the LFC files, pval files, W-score files all into parent folders with those titles. 

#Test this out for one qiime2 object. 
wscoretest <- read_csv("~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedWscores/ps1014hictrl-w_slice.csv")
wscoretest1 <- wscoretest[,c(1,3)]
colnames(wscoretest1)[2] <- "W"
lfctest <- read_csv("~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedLFCs/ps1014hictrl-lfc_slice.csv")
lfctest1 <- lfctest[,c(1,3)]
colnames(lfctest1)[2] <- "LFC"
pvaltest <- read_csv("~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedpvals/ps1014hictrl-p_val_slice.csv")
pvaltest1 <- pvaltest[,c(1,3)]
colnames(pvaltest1)[2] <- "pval"

ps1014_hictrl <- merge(pvaltest1, wscoretest1, by="id", all.x=TRUE)
ps1014_hictrl <- merge(ps1014_hictrl, lfctest1, by = "id", all.x=TRUE)

head(ps1014_hictrl)



#For all the Wscore files. 
# 1. Specify the path to your folder containing the CSV files
folder_path <- "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedWscores" # Replace with your actual path
# 2. Get a list of all file paths ending with ".tsv"
files_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
# 3. Read all files into a single data frame including adding a column to identify which file each row came from
new_col_name <- "ANCOMBC.W.score"
list_of_dfs <- lapply(files_list, function(file) {
  df <- read.csv(file) # Use read.table or another appropriate function for your file type
  names(df)[3] <- new_col_name # Rename the second column by index
  df$source <- basename(file) 
  return(df)
})
ancombcWscores <- do.call(rbind, list_of_dfs)
ancombcWscores$source <- sub("\\-[^-]+$", "", ancombcWscores$source) #This is a "greedy operator" that drops the last period and everything after it. s
ancombcWscores <- ancombcWscores[,c(1,3,4)]
ancombcWscores$bind <- paste(ancombcWscores$source, ancombcWscores$id, sep="_")
head(ancombcWscores)
#Write it out if you want to save the intermediate. 
write_xlsx(x = ancombcWscores, path = "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedWscores/ancombcWscores.xlsx") 
ancombcWscores_files <- data.frame(files_list)
write_xlsx(ancombcWscores_files, "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedWscores/ancombcWscores_files.xlsx")
#

#For all the LFC files. 
# 1. Specify the path to your folder containing the CSV files
folder_path <- "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedLFCs" # Replace with your actual path
# 2. Get a list of all file paths ending with ".tsv"
files_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
# 3. Read all files into a single data frame including adding a column to identify which file each row came from
new_col_name <- "ANCOMBC.LFC"
list_of_dfs <- lapply(files_list, function(file) {
  df <- read.csv(file) # Use read.table or another appropriate function for your file type
  names(df)[3] <- new_col_name # Rename the second column by index
  df$source <- basename(file) 
  return(df)
})
ancombcLFCs <- do.call(rbind, list_of_dfs)
ancombcLFCs$source <- sub("\\-[^-]+$", "", ancombcLFCs$source) #This is a "greedy operator" that drops the last period and everything after it. s
ancombcLFCs <- ancombcLFCs[,c(1,3,4)]
ancombcLFCs$bind <- paste(ancombcLFCs$source, ancombcLFCs$id, sep="_")
head(ancombcLFCs)
#Write it out if you want to save the intermediate. 
write_xlsx(x = ancombcLFCs, path = "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedLFCs/ancombcLFCs.xlsx") 
ancombcLFCs_files <- data.frame(files_list)
write_xlsx(ancombcLFCs_files, "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedLFCs/ancombcLFCs_files.xlsx")
#

#For all the p-value files. 
# 1. Specify the path to your folder containing the CSV files
folder_path <- "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedpvals" # Replace with your actual path
# 2. Get a list of all file paths ending with ".tsv"
files_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
# 3. Read all files into a single data frame including adding a column to identify which file each row came from
new_col_name <- "ANCOMBC.p.val"
list_of_dfs <- lapply(files_list, function(file) {
  df <- read.csv(file) # Use read.table or another appropriate function for your file type
  names(df)[3] <- new_col_name # Rename the second column by index
  df$source <- basename(file) 
  return(df)
})
ancombcpvalues <- do.call(rbind, list_of_dfs)
ancombcpvalues$source <- sub("\\-[^-]+$", "", ancombcpvalues$source) #This is a "greedy operator" that drops the last period and everything after it. s
ancombcpvalues <- ancombcpvalues[,c(1,3,4)]
ancombcpvalues$bind <- paste(ancombcpvalues$source, ancombcpvalues$id, sep="_")
head(ancombcpvalues)
#Write it out if you want to save the intermediate. 
write_xlsx(x = ancombcpvalues, path = "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedpvals/ancombcpvalues.xlsx") 
ancombcpvalues_files <- data.frame(files_list)
write_xlsx(ancombcpvalues_files, "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/amassedpvals/ancombcpvalues_files.xlsx")
#
#
#
#   ##### Renamelabels for comparison ID's and bind everything together #####
#Subset to only the columns we want to keep in the merged df. 
ancombcLFCs_sub <- ancombcLFCs[,c(2,4)]
ancombcWscores_sub <- ancombcWscores[,c(2,4)]

ancombcalltmp <- merge(ancombcpvalues, ancombcLFCs_sub, by = "bind", all = TRUE)
ancombcalltmp1 <- merge(ancombcalltmp, ancombcWscores_sub, by = "bind", all = TRUE)

#Replace the names of the files with the list file names.
ancombcalltmp2 <- ancombcalltmp1 %>% mutate(source = recode(source, 
                                                          "ps1014hitreat" = "HI_worm_1014",
                                                          "ps1014hictrl" = "HI_sham_1014",
                                                          "ps1014loctrl" = "LO_sham_1014",
                                                          "ps1014lotreat" ="LO_worm_1014",
                                                          "ps110hictrl"  ="HI_sham_110",
                                                          "ps110hitreat" ="HI_worm_110",
                                                          "ps110loctrl" ="LO_sham_110",
                                                          "ps110lotreat" ="LO_worm_110",
                                                          "ps1428hictrl" ="HI_sham_1428", 
                                                          "ps1428hitreat" ="HI_worm_1428", 
                                                          "ps1428loctrl" ="LO_sham_1428",
                                                          "ps1428lotreat" ="LO_worm_1428",
                                                          "ps2841hitreat" ="HI_worm_2841",
                                                          "ps2841loctrl" ="LO_sham_2841",
                                                          "ps2841lotreat" ="LO_worm_2841",
                                                          "pspre14hictrl" ="HI_sham_pre14",
                                                          "pspre14loctrl" ="LO_sham_pre14",
                                                          "pspre14lotreat" ="LO_worm_pre14",
                                                          "pspre1hictrl" ="HI_sham_pre1",
                                                          "pspre1hitreat" ="HI_worm_pre1",
                                                          "pspre1loctrl" ="LO_sham_pre1",
                                                          "pspre1lotreat" ="LO_worm_pre1",
                                                          "pspre41hictrl" ="HI_sham_pre41",
                                                          "pspre41hitreat" ="HI_worm_pre41",
                                                          "pspre41loctrl" ="LO_sham_pre41",
                                                          "pspre41lotreat" ="LO_worm_pre41"
                                                          ))
unique(ancombcalltmp2$source)
ancombcall <- separate(ancombcalltmp2, 
                        col=source, 
                        into=c("microb", "treatment", "comparison"),
                        sep="_")
head(ancombcall)

write.csv(ancombcall, "~/Documents/1U4U_16S/GG2_reclassification/qiime_diffabund/ancombc_outputs/exported_data/ancombc-all_concat.csv")
View(ancombcall)
#









#

#   ##### Merge the dataframes for DESeq2 and ANCOM #####
# For some reason the ANCOM results have truncated the ASV ID a character or two short. Cut
# all ASV ID's at a certain number of characters for both ANCOM and DESeq results so that everything merges well. 
ancombcall$shortASV <- strtrim(ancombcall$id, 22)
ResistDESeq$shortASV <- strtrim(rownames(ResistDESeq), 22)
ResilDESeq$shortASV <- strtrim(ResilDESeq$ASV, 22)

#Ancomdat prep
ancombcall$mergecol <- paste(ancombcall$microb, 
                              ancombcall$treatment, 
                              ancombcall$comparison, 
                              ancombcall$shortASV, 
                              sep="_")
#Trim to only the columns you want to combine with the DESeq results. 
ancombcmerge <- ancombcall[,c("ANCOMBC.LFC", "ANCOMBC.p.val", "ANCOMBC.W.score", "mergecol")]

#Create a new column to merge with the ANCOM data. 
ResistDESeq$mergecol <- paste(ResistDESeq$MT, ResistDESeq$comparison, ResistDESeq$shortASV, sep="_")
ResilDESeq$mergecol <- paste(ResilDESeq$MT, ResilDESeq$comparison, ResilDESeq$shortASV, sep="_")
#Make sure these mergecol's both agree. 
head(ResilDESeq$mergecol)
head(ancombcall$mergecol)
#
ancombcmergetest <- ancombcmerge[ancombcmerge$ANCOMBC.p.val <= 0.05,]
#Merge the dataframes by the custom merge column. Keep only columns that are represented by rows in both dfs. 
ResistDiffAbundallBC <- merge(ResistDESeq, ancombcmergetest, by = "mergecol", all = TRUE)
ResilDiffAbundallBC <- merge(ResilDESeq, ancombcmergetest, by = "mergecol", all = TRUE)
#
ResistDiffAbundBC <- merge(ResistDESeq, ancombcmergetest, by = "mergecol", all = FALSE)
ResilDiffAbundBC <- merge(ResilDESeq, ancombcmergetest, by = "mergecol", all = FALSE)
#

#Write those out to keep them for future. 
write.csv(ResistDiffAbundBC, "~/Documents/1U4U_16S/GG2_reclassification/DESeq2/Resistance_DESeqANCOMBC_concatenated.csv")
write.csv(ResilDiffAbundBC, "~/Documents/1U4U_16S/GG2_reclassification/DESeq2/Resilience_DESeqANCOMBC_concatenated.csv")
#
#

#     ##### UPDATED -- Horizontal barplots   #####
#Create a new column for naming colors and order MT columns. :
ResistDiffAbundBC$MTD <- paste(ResistDiffAbundBC$MT, ResistDiffAbundBC$comparison, sep = "_")
ResistDiffAbundBC$MTD <- factor(ResistDiffAbundBC$MTD, 
                              levels=c("LO_worm_110","LO_worm_2841", "HI_worm_110", 
                                       "HI_worm_1014", "HI_worm_1428","HI_worm_2841", 
                                       "LO_sham_1014", "HI_sham_110", "HI_sham_1014", "HI_sham_2841"))

ResistDiffAbundBC$MT <- factor(ResistDiffAbundBC$MT, levels=c("LO_worm","LO_sham", "HI_worm", "HI_sham"))

ResilDiffAbundBC$MTD <- paste(ResilDiffAbundBC$MT, ResilDiffAbundBC$comparison, sep = "_")
ResilDiffAbundBC$MTD <- factor(ResilDiffAbundBC$MTD, 
                             levels=c("LO_worm_pre14", "LO_worm_pre41", "HI_worm_pre41",
                                      "LO_sham_pre14", "LO_sham_pre41","HI_sham_pre14"))
ResilDiffAbundBC$MT <- factor(ResilDiffAbundBC$MT, levels=c("LO_worm","LO_sham", "HI_worm", "HI_sham"))
#Set up colors for 4-color calls based on MT. 
mycols <- c("LO_worm" = '#C66EF5' ,
            "LO_sham" = '#FFFFFF', 
            "HI_worm" = '#FDD964', 
            "HI_sham" = '#FFFFFF')

min(ResistDiffAbundBC$log2FoldChange)
max(ResistDiffAbundBC$log2FoldChange)

#        ##### ResistDiffAbundBC with all samples #####
min(ResistDiffAbundBC$ANCOMBC.LFC)
DFtoPlot <- ResistDiffAbundBC
#To make final figures, should probably subset so that faceted ASVs are all ordered correctly. 

#Plotting from ResistDiffAbundBC for worm-treated animals only. 
pDiffResist <- ggplot(DFtoPlot, aes(x=reorder(shortASV, -log2FoldChange), y =log2FoldChange)) +  
  #pDiffResist <- ggplot(DFtoPlot, aes(x=name1, y =log2FoldChange)) +
  geom_bar(stat = "identity", color= '#000000', aes(fill=MT)) +
  theme_bw() +
  scale_fill_manual(values = mycols) +
  #ylim(-6,8) +
  ggtitle("Resistance Timepoints Diff Abund Final") + 
  ylab("Log2 Fold Change") +
  xlab("ASV") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none") + 
  facet_wrap(~MTD, scales="free_y", space= "free_y", ncol=1) 
pDiffResistRotate <- pDiffResist  + coord_flip()
pDiffResistRotate
# values = setNames(c('red','green')  makes it so that positive is red. 

ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/ResistDiffAbundBC_Barplot.pdf",
       pDiffResistRotate, units="in", height=30, width=12)


#        ##### ResistDiffAbundBC with only worm-treated animals #####
min(ResistDiffAbundBC$ANCOMBC.LFC)
DFtoPlot <- ResistDiffAbundBC[ResistDiffAbundBC$MT %in% c("LO_worm", "HI_worm"),]
#To make final figures, should probably subset so that faceted ASVs are all ordered correctly. 

#Plotting from ResistDiffAbundBC for worm-treated animals only. 
pDiffResist <- ggplot(DFtoPlot, aes(x=reorder(name, -log2FoldChange), y =log2FoldChange)) +  
#pDiffResist <- ggplot(DFtoPlot, aes(x=name1, y =log2FoldChange)) +
    geom_bar(stat = "identity", color= '#000000', aes(fill=MT)) +
  theme_bw() +
  scale_fill_manual(values = mycols) +
  #ylim(-6,8) +
  ggtitle("Resistance Timepoints Diff Abund Final") + 
  ylab("Log2 Fold Change") +
  xlab("ASV") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none") + 
  facet_wrap(~MTD, scales="free_y", space= "free_y", ncol=1) 
pDiffResistRotate <- pDiffResist  + coord_flip()
pDiffResistRotate
# values = setNames(c('red','green')  makes it so that positive is red. 

ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/ResistDiffAbundBC_WormTreated_Barplot.pdf",
       pDiffResistRotate, units="in", height=30, width=12)

#        ##### ResistDiffAbundBC for 1 dpi to 10 dpi #####
DFtoPlot <- ResistDiffAbundBC[ResistDiffAbundBC$comparison =="110",]
#To make final figures, should probably subset so that faceted ASVs are all ordered correctly. 

#Plotting from ResistDiffAbundBC for worm-treated animals only. 
pDiffResist110 <- ggplot(DFtoPlot, aes(x=reorder(shortASV, -log2FoldChange), y =log2FoldChange)) +  
  geom_bar(stat = "identity", color= '#000000', aes(fill=MT)) +
  theme_bw() +
  scale_fill_manual(values = mycols) +
  #ylim(-6,8) +
  ggtitle("Diff Abund 1 dpi to 10 dpi") + 
  ylab("Log2 Fold Change") +
  xlab("ASV") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none") #+ facet_wrap(~MTD, scales="free_y", space= "free_y", ncol=1) 
pDiffResistRotate110 <- pDiffResist110  + coord_flip()
pDiffResistRotate110
# values = setNames(c('red','green')  makes it so that positive is red. 

ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/ResistDiffAbundBC_110comparison_notfaceted_Barplot.pdf",
       pDiffResistRotate110, units="in", height=25, width=10)


#        ##### Plotting from ResilDiffAbundBC #####
pDiffResil <- ggplot(ResilDiffAbundBC, aes(x=reorder(name, -log2FoldChange), y =log2FoldChange)) +  
  geom_bar(stat = "identity", color= '#000000', aes(fill=MT)) +
  theme_bw() +
  scale_fill_manual(values = mycols) +
  ylab("Log2 Fold Change") +
  xlab("ASV") +
  #ylim(-6,8) +
  ggtitle("Resilience Timepoints Diff Abund Final") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none") +
  facet_wrap(~MTD, scales="free_y", space= "free_y", ncol=1) 
pDiffResilRotate <- pDiffResil +  coord_flip()
pDiffResilRotate
#
ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/ResilDiffAbundBC_Barplot.pdf",
       pDiffResilRotate, units="in", height=13, width=10)

#


#   ##### Summarize data for numbers of differentially abundant ASVs per group, per time point #####
View(ResistDiffAbundBC)
unique(ResistDiffAbundBC$MTD)
#HI_sham_1014 HI_sham_110  HI_worm_1014 HI_worm_110  HI_worm_1428 HI_worm_2841 LO_sham_1014 LO_worm_110  LO_worm_2841

nrow(ResistDiffAbundBC[ResistDiffAbundBC$MTD == 'LO_worm_2841' & ResistDiffAbundBC$log2FoldChange > 0,])
nrow(ResistDiffAbundBC[ResistDiffAbundBC$MTD == 'LO_worm_2841' & ResistDiffAbundBC$log2FoldChange < 0,])

#LO_worm
##LO_worm_110: 50 up, 14 down
##LO_worm_2841: 2 up, 0 down

#LO_sham
##LO_sham_1014: 0 up, 1 down

#HI_worm
##HI_worm_110: 16 up, 7 down
##HI_worm_1014: 2 up, 23 down
##HI_worm_1428: 0 up, 6 down
##HI_worm_2841: 3 up, 0 down

#HI_sham
##HI_sham_110: 12 up, 16 down
##HI_sham_1014:0 up, 9 down
##

#
#   ##### Horizontal barplot like Lefse #####
#Create a new column for naming colors and order MT columns. :
ResistDiffAbund$MTD <- paste(ResistDiffAbund$MT, ResistDiffAbund$comparison, sep = "_")
ResistDiffAbund$MTD <- factor(ResistDiffAbund$MTD, 
                              levels=c("LO_worm_110","LO_worm_2841", "HI_worm_110", 
                                       "HI_worm_1014", "HI_worm_1428","LO_sham_1014",
                                       "HI_sham_2841"))
ResistDiffAbund$MT <- factor(ResistDiffAbund$MT, levels=c("LO_worm","LO_sham", "HI_worm", "HI_sham"))

ResilDiffAbund$MTD <- paste(ResilDiffAbund$MT, ResilDiffAbund$comparison, sep = "_")
ResilDiffAbund$MTD <- factor(ResilDiffAbund$MTD, 
                             levels=c("LO_worm_pre14", "LO_worm_pre41", "HI_worm_pre41",
                                      "LO_sham_pre14", "LO_sham_pre41","HI_sham_pre14"))
ResilDiffAbund$MT <- factor(ResilDiffAbund$MT, levels=c("LO_worm","LO_sham", "HI_worm", "HI_sham"))
#Set up colors for 4-color calls based on MT. 
mycols <- c('#C66EF5' ,'#FFFFFF', '#FDD964', '#FFFFFF')

min(ResistDiffAbund$log2FoldChange)
max(ResistDiffAbund$log2FoldChange)

#   ##### Make a dot plot with taxonomy.  #####
#Now to look at the differentially enriched taxa in ggplot
#library("ggplot2")

x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#Plot it at family level
X <- ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Family)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ggtitle("DESeq2 pspre14lotreat_fam ~ day w/lfcShrink")
X

#Double check which way is increase and decrease log fold change. 
#I think that when comparing D14_vs_D.5, enrichment refers to 1st -> 2nd, so down means
# fewer reads in 2nd group (in this case, Dpre).



#   ##### Plotting counts of each microbe #####
#See what the counts look like for the different metabolites.
#Either using the defaults for axis scales
plotCounts(dsps, gene=which.min(resLFC$padj), intgroup="day")

#Or to customize the plot in ggplot. 
d <- plotCounts(dds, gene="C733...4.95", intgroup="label", 
                returnData=TRUE) #returns a df that can be used with ggplot.
ggplot(d, aes(x=label, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0))

#which.min(res$padj)
#scale_y_log10(breaks=c(25,100,400))




