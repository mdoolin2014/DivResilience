#Loading in microbiome data and subsetting functions

# Updated from October 2021, in March 2022 and again several times with updated data outputs from QIIME2 until July 2023. 
# Adapted from previous versions from Dylan Klure (2019) that started from raw sequences and from the tutorials: 
# https://benjjneb.github.io/dada2/tutorial.html
# https://joey711.github.io/phyloseq/
# Make sure all packages and R is up-to-date first! Comments are my own.

#Importing 1U4U qiime2 outputs using qiime2R and working from there. Imports created 
# using qiime2-2022.11 on the CHPC (previously using 2021.4). SLURM script 1U4U_reclassify_16S_slurm.sh 
# Bsed deblur denoising, classified with V3V4 region of Greengenes2 database with backbone 
# made by Zac Stephens for better resolution and more accurate UniFrac results. 
#This whole dataset includes samples from LO bubble presamples, HI and LO experimental
# samples, and run controls. Will need to filter those down.

setwd("~/Documents/1U4U_16S/GG2_reclassification/")

locols <- c('#FDF8FF' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
hicols <- c("#FFFCF3" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")

library(ggplot2)
library(ape)
library(gridExtra)
library(SummarizedExperiment)
library(dada2) 
library(S4Vectors)
library(stats4)
library(IRanges)
library(XVector)
library(RSQLite)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(vegan)
library(reshape2)
library(phyloseq)
library(Rcpp)

#install all dependencies, which are a lot since I just updated my R. 
library(qiime2R) #package version 0.99.6

#Make sure you're in the right working directory to pull qiime outputs. 
setwd("~/Documents/1U4U_16S/GG2_reclassification")
#qza_to_phyloseq(features, tree, taxonomy, metadata, tmp)

#Loading in my physeq object. I've already filtered out chloroplast and mitochondria
# reads, as well as singletons and doubletons, in qiime. Can go ahead and make ps1.
ps1 <- qza_to_phyloseq(features="final-table.qza", taxonomy="taxonomy_gg2.qza", tree ="tree_insertion_root.qza",
                       metadata="~/Documents/1U4U_16S/Rerun_metadata/MapBySampleID_1U4U_w18684.tsv")

saveRDS(ps1, "ps1RDS.rds")


#Subsetting my original phyloseq object to have only samples and controls
ps2 <- subset_samples(ps1, sample_type=="sample")
ps2 <- prune_taxa(taxa_sums(ps2) >= 1, ps2) #Take the missing data out for things that matter

#psexp is unrarefied datasets with all experimental days. 
psexptmp <- subset_samples(ps2, !day %in% c("none", "pre") & mouse!="18_RLP")
psexptmp <- prune_taxa(taxa_sums(psexptmp)>=1, psexptmp) #This includes 345 samples
psexptmpno31 <- subset_samples(psexptmp, day != "D31")
psexptmpno31 <- prune_taxa(taxa_sums(psexptmpno31)>=1, psexptmpno31)
mean(sample_sums(psexptmpno31))
##Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##801    7485   18371   21010   29530   97519 , sd = 15976.32
samlist <- sample_sums(psexptmp) > 3000
psexp <- prune_samples(samlist, psexptmp)
psexp #ok now we have the same samples in the unrarefied data as psexp3000
psexp@sam_data[["day_num"]] <- as.numeric(gsub("D", "", psexp@sam_data[["day"]]))
#324 samples left after this trim and removing 18_RLP because it's weird. 
samlist <- sample_sums(psexptmpno31) > 3000
psexpno31 <- prune_samples(samlist, psexptmpno31)
#

psL3 <- subset_samples(ps1, sample_group == "L3")
psL3 <- prune_taxa(taxa_sums(psL3) >=1, psL3)


#Pull metadata
library(devtools)
#devtools::install_github('twbattaglia/btools')
library(btools) #version 0.0.1, installed on 14Oct2021
obj <- psexp3000no31
#Make dataframe from sample data. 
dftmp <- data.frame(estimate_richness(obj, measures = "Observed"),
                    sum = sample_sums(obj), 
                    estimate_richness(obj, measures="Shannon"),
                    btools::estimate_pd(obj),    #This adds Faith's phylo div and spec richness
                    "row"=rownames(sample_data(obj)),
                    evenness=(estimate_richness(obj, measures="Shannon")/log(estimate_richness(obj, measures = "Observed"))),
                    "cage_day"=paste0(sample_data(obj)$cage,sample_data(obj)$day))  
#Don't forget to update final df name as well. 
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
dfexp3000no31<- cbind(sample_data(obj),dftmp)
head(dfexp3000no31)

#


# Going from R back into qiime #
#Start with a phyloseq object.
library(biomformat) #version 1.28.0
#Creating a phyloseq object that only contains the samples that are also in the 
# meatabolomics dataset
samps <- data.frame(read.csv("metab_sample_list.csv")) #Get the sample names. 
#Create the taxa table for use as FeatureTable[Sequence] (I think)
tax <- as(tax_table(psexp3000no31), "matrix")
tax_cols <- colnames(tax)
tax <- as.data.frame(tax)
tax$taxonomy <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co] <- NULL
write.table(tax, "~/Documents/psexp3000no31_tax.tsv", quote=FALSE, col.names=FALSE, sep="\t")

#Create your OTU table for the biom file, as FeatureTable[Frequency]
otu <- as(otu_table(psdecont),"matrix") #taxa_are_rows, so don't need to transpose.
otu_biom <- make_biom(data=otu)
write_biom(otu_biom, "~/Documents/psexp3000no31.biom")

#Could also export this as a text file, but biom file from phyloseq object
# seems easier.

#Now write out your metadata for this phyloseq object.
write.table(data.frame(sample_data(psexp3000no31)), "psdecon.tsv", sep="\t",
            row.names=TRUE, col.names=TRUE, quote=FALSE)

#