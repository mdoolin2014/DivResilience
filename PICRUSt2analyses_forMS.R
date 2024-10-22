#PICRUSt2 outputs analyses for "Low diversity leads to instability and resilience following perturbation"
# submitted to The ISME Journal. 
#Full PICRUSt2 pipeline was run in terminal, then output count table and KO metadata were imported  here for analyses.
#September 2024, M.L. Doolin


##### Load libraries, set up workspace, hex codes, etc. #####

#Libraries will maybe already be loaded, but just in case
library(easypackages) #Package to be able to load multiple libraries in one command
libraries("ggplot2", "ggpubr", "gridExtra", "ggbreak", "ggbeeswarm", "ggstatsplot", "ggpicrust2", "ggfortify", "scales") #plotting packages
libraries("lefser", "ALDEx2", "edgeR", "Maaslin2", "ComplexHeatmap") #special differntial abundance stuff
library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(readr)
library(reshape2)
library(tidyr)
library(funrar)
library(vegan)
library('stringr')
#library(ecodist) #might not need this anymore, can't remember why I had it in.
library(ape)

tax_table <- phyloseq::tax_table

locols <- c('#FDF8FF' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
hicols <- c("#FFFCF3" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")


loimpcols <- c('#FDF8FF' , '#A11FE5', '#300049')
hiimpcols <- c("#FFFCF3" , "#DCAE1C" , "#644C00")

loimpellipse <- c('#c2b6cf' , '#A11FE5', '#300049')
hiimpellipse <- c("#c9c4b3" , "#DCAE1C" , "#644C00")

#If you need to make sure the boxplot for D-5 can be seen, change the first color to black. 
otherlo <-  c('#999999' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
otherhi <- c("#999999" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")


#Colors to use based on microbiome and timepoint, with variable "MT"
allcols <- c("#FDF8FF",  '#d6a3f0' , '#c275eb' , '#A11FE5' , '#7E05BE' , '#5C008D', '#300049', 
             "#faf1d7", "#f5e1a6" , "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")
#Colors to use to define all control samples as grey, with variable being "MDT". 
allcols2 <- c("#e0c0f0",  '#d6a3f0' , '#c275eb' , '#A11FE5' , '#7E05BE' , '#5C008D', '#300049',
              "#faf1d7", "#f5e1a6" , "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00", 
              "#6c696e", "#6c696e", "#6c696e", "#6c696e",  "#6c696e", "#6c696e", "#6c696e",
              "#6c696e", "#6c696e", "#6c696e", "#6c696e",  "#6c696e", "#6c696e", "#6c696e")
allcols1 <- c("#e0c0f0",  '#d6a3f0' , '#c275eb' , '#A11FE5' , '#7E05BE' , '#5C008D', '#300049',
              "#faf1d7", "#f5e1a6" , "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00", 
              "#FFFFFF", "#FFFFFF", "#FFFFFF",  "#FFFFFF", "#FFFFFF", "#FFFFFF","#FFFFFF", "#FFFFFF", 
              "#FFFFFF", "#FFFFFF", "#FFFFFF",  "#FFFFFF", "#FFFFFF", "#FFFFFF")
#





###### DID USE FOR 1U4U Create a KO phyloseq object to do phyloseq_to_deseq #####

#Use the unstrat dataset for this, without descriptions. This is what I'll convert to 
# read in as an OTU table for a phyloseq object. Plus, metadata
setwd("/Users/mdoolin")

#Pull in metadata.
META <- data.frame(read_tsv("/Users/mdoolin/psexp_metadat.tsv"))
rownames(META) <- META$sample.id
META <- as.data.frame(META)
META <- sample_data(META)

#Create OTU table from picrust output.
#I rounded these for diff abundance just to see if that would help. Could try not 
# to do this. 
OTU <- data.frame(read_tsv("pred_metagenome_unstrat.tsv"))
rownames(OTU) <- OTU$KEGG #This will work if your first column is named "KEGG"
OTU <- as.matrix(OTU[,-1]) #remove the label column when you make it a matrix. Only numeric left. 
OTU <- round(OTU)  #round counts to whole integers for subsequent processing. 
colnames(OTU) <- gsub('X', '', as.character(colnames(OTU)))
OTU <- otu_table(OTU, taxa_are_rows=TRUE)


#Make a copy of the KO description file and create a taxonomy table from it. 
TAX <- data.frame(read_tsv("/Users/mdoolin/KEGG_functions.tsv"))
nrow(TAX)
TAX$B <- str_replace_all(TAX$B,':',' -')
TAX$C <- str_replace_all(TAX$C,':',' -')
TAX1 <- data.frame('KEGG' = TAX$KEGG, 
                   'func' = paste(TAX$A, TAX$B, TAX$C, sep=": "))
TAX1$KEGG <- make.unique(TAX1$KEGG) #Need to make unique names for each ko rowname for casting (opposite of melt)
TAX2 <- TAX1 %>% separate_wider_delim(KEGG, delim= ".", names=c('ko', 'replicate'), too_few = "align_start") #separate the column based on the .
#use "too_few="align_start"" bc some of the numbers (the first in each series) will not have a .#
TAX2 <- data.frame(TAX2) #using tidyr function made this into a tibble. Convert back to df. 
TAX2[is.na(TAX2)] <- 0 #convert the NAs in the replicate column to 0 for a unique replicate name. 
TAX3 <- TAX2 %>% pivot_wider(names_from=replicate, values_from=func)#This is your complete DF. Let's just make 
# a df with the first recognized functions for simplicity. 
TAX4 <- data.frame('KEGG' = TAX3$ko, 
                   'func' = TAX3$`0`) 
TAX5 <- TAX4 %>% separate_wider_delim(func, 
                                      delim= ":", 
                                      names=c('General', 'Specific', 'name'), 
                                      too_many="merge") #separate the column based on the :
TAX5 <- data.frame(TAX5)
rownames(TAX5) <- TAX5$KEGG #hmm, this is an issue because there are a bunch of duplicates since 
# many of these orthologs are involved in many processes. 
TAX5 <- TAX5[,-1]
TAX5 <- as.matrix(TAX5)
TAX <- tax_table(TAX5)


#   ##### Put together the pieces to make a phyloseq object. #####

#Create your phyloseq object from the GAT DAMN pieces. 
psko <- phyloseq(OTU, TAX, sample_data=META)
psko

#


sample_data(psko)$MDT <- paste(sample_data(pstrans1)$div_day, sample_data(pstrans1)$treatment, sep="_")
sample_data(psko)$microb <- factor(sample_data(psko)$microb, levels=c("LO", "HI"))
sample_data(psko)$treatment <- factor(sample_data(psko)$treatment, levels=c("control", "worm"))
sample_data(psko)$div_day <- factor(sample_data(psko)$div_day, levels=c('LO_D-5', 'LO_D1', 'LO_D10', 'LO_D14', 'LO_D28', 
                                                                        'LO_D31', 'LO_D41', 'HI_D-5', 'HI_D1', 'HI_D10', 
                                                                        'HI_D14', 'HI_D28', 'HI_D31', 'HI_D41'))
sample_data(psko)$MDT <- factor(sample_data(psko)$MDT, 
                                levels=c('LO_D-5_worm', 'LO_D1_worm', 'LO_D10_worm', 'LO_D14_worm', 'LO_D28_worm','LO_D31_worm', 'LO_D41_worm',
                                         'HI_D-5_worm', 'HI_D1_worm', 'HI_D10_worm', 'HI_D14_worm', 'HI_D28_worm', 'HI_D31_worm', 'HI_D41_worm', 
                                         'LO_D-5_control', 'LO_D1_control', 'LO_D10_control', 'LO_D14_control', 'LO_D28_control', 'LO_D31_control', 'LO_D41_control', 
                                         'HI_D-5_control', 'HI_D1_control', 'HI_D10_control', 'HI_D14_control', 'HI_D28_control', 'HI_D31_control', 'HI_D41_control'))

sample_data(psko) %>% rename('ASVs' = 'Observed')
#

#   ##### Subset your phyloseq object into useful subsets. #####

#For pairwise comparisons

#psprelotreat
pstmp <- subset_samples(psko, day %in% c("D-5", "D41") & 
                          microb == "LO" & 
                          treatment == "worm")
pspre41lotreat <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4977 functions in 55 samples

#pspre41hitreat
pstmp <- subset_samples(psko, day %in% c("D-5", "D41") & 
                          microb == "HI" & 
                          treatment == "worm")
pspre41hitreat <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#3587 taxa in 26 samples

#pspre41loctrl
pstmp <- subset_samples(psko, day %in% c("D-5", "D41") & 
                          microb == "LO" & 
                          treatment == "control")
pspre41loctrl <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4217 taxa in 12 samples

#pspre41hictrl
pstmp <- subset_samples(psko, day %in% c("D-5", "D41") & 
                          microb == "HI" & 
                          treatment == "control")
pspre41hictrl <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#3496 taxa in 8 samples

#Only one day and one treatment group. 
#pspretreat
pstmp <- subset_samples(psko, day %in% c("D-5") & treatment == "worm")
pspretreat <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#5038 functions in 50 samples

#psprectrl
pstmp <- subset_samples(psko, day %in% c("D-5") & treatment == "control")
psprectrl<- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4411 taxa in 10 samples

#ps41treat
pstmp <- subset_samples(psko, day %in% c("D41") & treatment == "worm")
ps41treat<- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4634 taxa in 31 samples

#ps41ctrl
pstmp <- subset_samples(psko, day %in% c("D41") & treatment == "control")
ps41ctrl<- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4409 taxa in 10 samples

#Only one day and one microbiome type
#psprelo
pstmp <- subset_samples(psko, day %in% c("D-5") & microb == "LO")
psprelo <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4883 taxa in 38 samples

#psprehi
pstmp <- subset_samples(psko, day %in% c("D-5") & microb == "HI")
psprehi <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#3525 taxa in 22 samples

#ps41lo
pstmp <- subset_samples(psko, day %in% c('D41') & microb == "LO")
ps41lo <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#4471 taxa in 29 samples

#ps41hi
pstmp <- subset_samples(psko, day %in% c('D41') & microb == "HI")
ps41hi <- prune_taxa(taxa_sums(pstmp) > 1, pstmp)
#3554 taxa in 12 samples
#



##### Run deseq on this picrust phyloseq object. #####
library(DESeq2)
library(microbiome)
library(ashr)
library(ape)

#Tutorial on DESeq2 from phyloseq.
#https://joey711.github.io/phyloseq-extensions/DESeq2.html

#   ##### Run DESeq2 at the level of pathway. #####
obj <- pspre41hitreat

#sample_data(obj)$day <- as.factor(sample_data(obj)$day)
#Using phyloseq_to_deseq2 to do the hard work and transform data into a form
# that is able to be analyzed by this program, which was made for RNA seq data
dsps = phyloseq_to_deseq2(obj, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)

summary(res)

#Need to run LFCShrink
resultsNames(dsps)   # to get coef name
#names are: 'microb_LO_vs_HI', 'treatment_worm_vs_control', 'day_D41_vs_D.5'
#positive is enriched in the group that's listed first
resLFC <- lfcShrink(dsps, coef="day_D41_vs_D.5", type="ashr") #change coef
summary(resLFC)


resOrdered <- resLFC[order(resLFC$log2FoldChange),]
resNames <- rownames(resOrdered)
tmp <- data.frame(tax_table(psko))

#All outcomes, for volcano plots. 
deseq.all=cbind(data.frame(resOrdered), tmp[rownames(resOrdered), ])


#Now to only grab the significant results. 
a = 0.05
b = 1
sigtab = resOrdered[which(resOrdered$padj < a), ]
sigtab = sigtab[abs(sigtab$log2FoldChange) >= b, ]
sigtab = cbind(data.frame(sigtab), tmp[rownames(sigtab), ])
sigtab$final <- paste(sigtab$name, sep="  ", rownames(sigtab))
dim(sigtab)
View(sigtab)

#note that there were no sig diff for psprelo or psprehi by treatment, nor was there any
# difference for pspre41hitreat and pspre41hictrl by day. 
tab8 <- cbind(sigtab, "comparison" = "ps41hi") # 96 pass cutoffs
tab7 <- cbind(sigtab, "comparison" = "ps41lo") # 369 pass cutoffs
#tab6 <- cbind(sigtab, "comparison" = "ps41treat") # 2092 pass cutoffs
#tab5 <- cbind(sigtab, "comparison" = "ps41ctrl") # 1117 pass cutoffs
#tab4 <- cbind(sigtab, "comparison" = "psprectrl") #2154 pass cutoffs
#tab3 <- cbind(sigtab, "comparison" = "pspretreat") #2696 pass cutoff
tab2 <- cbind(sigtab, "comparison" = "pspre41loctrl") #318 pass cutoffs
tab1 <- cbind(sigtab, "comparison" = "pspre41lotreat") #2317 pass cutoffs

#sigtabfull <- rbind(tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8)
#write.csv(sigtabfull, "allsignif_lowestRank_DESeq2.csv")


#Now to look at the differentially enriched taxa in ggplot
#library("ggplot2")

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$General, function(x) max(x))
x = sort(x, TRUE)
sigtab$General = factor(as.character(sigtab$General), levels=names(x))



#Plot it at family level
X <- ggplot(sigtab, aes(x=General, y=log2FoldChange, color=General)) + geom_point(size=6) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ggtitle("DESeq2  ~ Day w/lfcShrink")
X



#Make bar graphs to summarize the relative amounts of each function that were 
df1 <- tab2[tab2$log2FoldChange>0,]
df2 <- tab2[tab2$log2FoldChange<0,]

path_counts_pos <- df1 %>%
  group_by(General) %>%
  summarise(n = n())

path_counts_neg <- df2 %>%
  group_by(General) %>%
  summarise(n = n())

#How many total pathways were differentially abundant?
sum(path_counts_pos$n)

ggp1 <- ggplot(path_counts_pos, aes(x=General, y=n)) + 
  geom_bar(aes(x=General, y=n), stat="identity") +
  rotate_x_text(40) +
  xlab("High-Level KO function") +
  ylab('Number of differential KO counts') +
  ggtitle('Enriched at D41 in PhyLo11B sham-treated')


#
#   ##### Combine functions at a higher level to look at broad-scale changes. #####

ps1 <- tax_glom(pspre41lotreat, taxrank=rank_names(pspre41lotreat)[2])
ps1

sample_data(ps1)$day <- as.factor(sample_data(ps1)$day)

#Using phyloseq_to_deseq2 to do the hard work and transform data into a form
# that is able to be analyzed by this program, which was made for RNA seq data
dsps = phyloseq_to_deseq2(ps1, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
summary(res)
#shrink the LFC
resultsNames(dsps)   # to get coef name
resLFC <- lfcShrink(dsps, coef="day_D41_vs_D.5", type="ashr") #change coef
# can use apeglm or ashr as type. Note that they give different LFC values
# that are several orders of magnitude different, even if total ASVs are the same.
summary(resLFC)

##LFCShrink results for pspre41lotreat by day


resOrdered <- resLFC[order(resLFC$log2FoldChange),]
resNames <- rownames(resOrdered)
tmp <- data.frame(tax_table(ps1))

#Now to only grab the significant results. 
a = 0.05
b = 1
sigtab = resOrdered[which(resOrdered$padj < a), ]
sigtab = sigtab[abs(sigtab$log2FoldChange) >= b, ]
sigtab = cbind(data.frame(sigtab), tmp[rownames(sigtab), ])
sigtab$final <- paste(sigtab$name, sep="  ", rownames(sigtab))
dim(sigtab)
View(sigtab)


#tab7 <- cbind(sigtab, "comparison" = "ps41lo")
#tab6 <- cbind(sigtab, "comparison" = "ps41treat")
#tab5 <- cbind(sigtab, "comparison" = "ps41ctrl")
#tab4 <- cbind(sigtab, "comparison" = "psprectrl")
#tab3 <- cbind(sigtab, "comparison" = "pspretreat")
#tab2 <- cbind(sigtab, "comparison" = "pspre41loctrl")
#tab1 <- cbind(sigtab, "comparison" = "pspre41lotreat")

#sigtabfull <- rbind(tab1, tab2, tab3, tab4, tab5, tab6, tab7)
#write.csv(sigtabfull, "allsignif_midlevelRank_DESeq2.csv")


#      ##### Dot plot of the DESeq2 results #####

#Now to look at the differentially enriched taxa in ggplot
#library("ggplot2")

sigtab <- read.csv("allsignif_midlevelRank_DESeq2.csv")

tab1 <- subset(sigtab, comparison %in% c("pspre41lotreat", "pspre41loctrl"))
tab1$factor <- "day"

tab2 <- subset(sigtab, comparison %in% c("pspretreat", "psprectrl"))
tab2$factor <- "microbpre"

tab3 <- subset(sigtab, comparison %in% c("ps41lo"))
tab3$factor <- "treatment"

tab4 <- subset(sigtab, comparison %in% c("ps41ctrl", "ps41treat"))
tab4$factor <- 'microbpost'

fulltab <- rbind(tab1, tab2, tab3, tab4)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set2", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Order the KO functions by the highest level. 
x = tapply(fulltab$log2FoldChange, fulltab$General, function(x) max(x))
x = sort(x, TRUE)
fulltab$General = factor(as.character(fulltab$General), levels=names(x))


#Plot it at family level
X <- ggplot(fulltab, aes(x=Specific, y=log2FoldChange, color=comparison)) + geom_point(size=6) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  ggtitle("DESeq2 at High-level Function Rank w/lfcShrink") +
  facet_wrap(~ factor, scales= "free_x", ncol=4)
X


ggsave("DESeq2_allresultsfaceted_MidlevelRankGlom.pdf", X, units="in", height=6, width=13)


#

##### Beta diversity analyses for ISME #####
#   ##### Betadisper #####
library(vegan) #version 2.5-7

#Before deciding which beta div analysis, run a dispersion test to see if appropriate for parametrics
x <- phyloseq::distance(psko, 'bray')
y <- as(sample_data(psko), "data.frame")
groups <- y[["microb"]]
mod <- betadisper(x, groups)
anova(mod)
plot(mod)  #To see what the points around the centroid look like
plot(TukeyHSD(mod))
permutest(mod)


#   ##### Create distance matrices (pstrans iterations) and PCoAs for psexp#####
pstrans1 <- transform_sample_counts(psko, function(otu) otu/sum(otu)) 
#Create ordinations.
ps.ordi <- ordinate(pstrans1, "PCoA", "bray")
ps.ordi2 <- ordinate(pstrans1, "PCoA", "jaccard")
plot_scree(ps.ordi2) #for PCoA, to plot eigenvalues

#


brayp <- plot_ordination(pstrans1, ps.ordi) +
  geom_point(aes(colour=microb, fill=MDT, shape=microb), size=3) +
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values=allcols1) +
  #scale_colour_manual(values=c("#62148c","#62148c",'#b58e14', "#b58e14")) +
  scale_colour_manual(values=c('#62148c', '#b58e14')) +
  #scale_colour_manual(values=c('#FFFFFF', '#62148c', '#FFFFFF', '#b58e14')) +
  # ylim(-0.45, 0.4)+
  # xlim(-0.63, 0.4)+
  stat_ellipse(aes(color= microb), type="norm", level=0.95)+
  theme_bw() +
  ggtitle("PICRUSt2 Bray PCOA") +
  theme(plot.title=element_text(size=16, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        axis.text=element_text(size=15), 
        legend.position="none")
brayp

ggsave("plots/PICRUSt2_allexpsamps_BrayPCOA.pdf", brayp, units="in", height=6, width=7)



#   ##### Create distance matrices and PCOAs for subsetted data. #####

#For PhyLo11B animals
pstrans1 <- transform_sample_counts(pspre41lotreat, function(otu) otu/sum(otu)) 
#Create ordinations.
ps.ordi1 <- ordinate(pstrans1, "PCoA", "bray")
ps.ordi2 <- ordinate(pstrans1, "PCoA", "jaccard")
plot_scree(ps.ordi2) #for PCoA, to plot eigenvalues


brayplo <- plot_ordination(pstrans1, ps.ordi1) +
  geom_point(aes(fill=day, shape=microb), size=3) +
  scale_shape_manual(values=c(21)) +
  scale_fill_manual(values=c('#e0c0f0', '#300049', '#e0c0f0', '#300049')) +
  stat_ellipse(aes(color= day), type="norm", level=0.95)+
  scale_colour_manual(values=c('#e0c0f0', '#300049')) +
  theme_bw() +
  ggtitle("PICRUSt2 Bray PCOA, pspre41lotreat") +
  theme(plot.title=element_text(size=16, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        axis.text=element_text(size=15))
brayplo



ggsave("plots/PICRUSt2_pspre41lotreat_BrayPCOA.pdf", brayplo, units="in", height=6, width=7)
#



#   ##### Run the PERMANOVAs #####

metadata <- as(sample_data(ps41hi), "data.frame")  
perms <- with(metadata, how(nperm = 5000)) #Note increasing number of permutations makes results more consistent.
#mod <- adonis2(phyloseq::distance(psexp3000lo, method="jaccard") ~  microb + day_num + microb*day_num, data=metadata) #Only one treatment
mod <- adonis2(phyloseq::distance(ps41hi, method="bray") ~ treatment, data=metadata) #Only one treatment
mod


#
##### Brief look at richness per sample. #####

obj <- pspre41hictrl

dat <- cbind(data.frame(sample_data(obj)), estimate_richness(obj, measures=c("Observed")))
names(dat)[44] <- "pathways" #renaming the picrust richness bc there is already ASV richness. 

summary(dat$pathways[dat$day=="D-5"])
sd(dat$pathways[dat$day=="D-5"])









