#Script for all the figures of the 1U4U paper. 


#####################################  PREPARING DATA  ###########################

############# Setting up workspace, loading in phyloseq object, subsetting objects #################

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

remotes::install_github("jbisanz/qiime2R")
library(qiime2R) #package version 0.99.6

setwd("~/Documents/1U4U_16S/GG2_reclassification")

#Loading in my physeq object from QIIME2. I've already filtered out chloroplast and mitochondria
# reads, as well as singletons and doubletons, in qiime. Can go ahead and make ps1.
ps1 <- qza_to_phyloseq(features="final-table.qza", taxonomy="taxonomy_gg2.qza", tree ="tree_insertion_root.qza",
                       metadata="~/Documents/1U4U_16S/Rerun_metadata/MapBySampleID_1U4U_w18684.tsv")

ps1@sam_data[["div_day"]] <- as.factor(paste(sample_data(ps1)$microb, sep="_", sample_data(ps1)$day))
ps1@sam_data[["TD"]] <- as.factor(paste(sample_data(ps1)$treatment, sep="_", sample_data(ps1)$day))
ps1@sam_data[["microb"]] <- factor(ps1@sam_data[["microb"]], levels=c("LO","HI"))
ps1@sam_data[["sam_day"]] <- as.factor(paste(sample_data(ps1)$treatment, sep="_", sample_data(ps1)$day))

View(sample_data(ps1))

ps2 <- subset_samples(ps1, sample_type=="sample")
ps2 <- prune_taxa(taxa_sums(ps2) >= 1, ps2) 

#Example of further subsetting based on sample_data. Then have subsetted into many different objects. 
psexp <- subset_samples(ps2, !day %in% c("none", "pre") & mouse!="18_RLP")
psexp <- prune_taxa(taxa_sums(psexp)>=1, psexp)
psexp
psexp@sam_data[["day_num"]] <- as.numeric(gsub("D", "", psexp@sam_data[["day"]]))


#Exporting as Taxa and OTU tables. 
write.csv(as.table(tax_table(psexp3000)), "psexp3000_TaxaTabGTDB.csv")
write.csv(as.table(otu_table(psexp3000)), "psexp3000_OTBTableGTSB.csv")
#Or concatenate them and export as single file. 
OTU <- data.frame(otu_table(psbubble))
OTU$sum <- rowSums(OTU)
otu_tax <- cbind(OTU, tax_table(psbubble), 
                 "ASV"=rownames(OTU))
#View(otu_tax)

write.csv(otu_tax, "psbubble_otutax_bind_sums.csv")


#Rarefy to 3000 reads. 
psexp3000 = rarefy_even_depth(psexp, rngseed = 52, sample.size=3000)
psexp3000 = subset_samples(psexp3000, mouse !="18_RLP") #Remove a mouse that looks really weird compared to everyone else. 
psexp3000 <- prune_taxa(taxa_sums(psexp3000)>= 1, psexp3000)

#Example of calling taxa names from a phyloseq object. 
nrow(table(phyloseq::tax_table(psexp3000)[,"Family"])) 


# Doing some QC on read depths #

tmpdf <- as.data.frame(sample_data(psexp))
df <- cbind(tmpdf, 
            "filtered"=sample_sums(psexp), 
            estimate_richness(psexp, measures = "Observed"))

boxplot(filtered ~ TD, data=df)
a <- aov(filtered~lane, data=df)
summary(a)
TukeyHSD(a)
plot(TukeyHSD(a))
#Original lane 2 has significantly higher read count after filtering, than
# Lane 1 or the rerun samples. 

############# Looking at your phyloseq objects, read counts, inspecting data ################


#What about for the rarefied dataset
tmpdf <- as.data.frame(sample_data(psexp))
df <- cbind(tmpdf, 
            "filtered"=sample_sums(psexp),
            estimate_richness(psexp, measures = "Observed"))

boxplot(Observed ~ lane, data=df)
a <- aov(Observed~lane, data=df)
summary(a)
TukeyHSD(a)
plot(TukeyHSD(a))


daysums <- data.frame("sample"=rownames(sample_data(ps2)), 
                      "mouse"=sample_data(ps2)$mouse,
                      "day"=sample_data(ps2)$day, 
                      "reads"=sample_sums(ps2))

lowsamples <- daysums[daysums$reads < 5000,]
View(lowsamples)


############# Pulling metadata from your phyloseq objects ##################
library(btools) #version 0.0.1, installed on 14Oct2021

obj <- psbubble

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
dfbubble<- cbind(sample_data(obj),dftmp)
head(dfexp)

#Subset df with metadata. 
dfexp3000treat <- dfexp3000[dfexp3000$treatment=="worm",]




#


############# Going from R back into QIIME2 ###############

#Start with a phyloseq object.

library(biomformat) #version 1.28.0

samps <- data.frame(read.csv("metab_sample_list.csv")) #Get the sample names that you want.  

#Create the taxa table for use as FeatureTable[Sequence] (I think)
tax <- as(tax_table(psexp3000), "matrix")
tax_cols <- colnames(tax)
tax <- as.data.frame(tax)
tax$taxonomy <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co] <- NULL
write.table(tax, "biom_inputs/psexp3000_tax.tsv", quote=FALSE, col.names=FALSE, sep="\t")

#Create your OTU table for the biom file, as FeatureTable[Frequency]
otu <- as(otu_table(psexp3000),"matrix") #taxa_are_rows, so don't need to transpose.
otu_biom <- make_biom(data=otu)
write_biom(otu_biom, "biom_inputs/psexp3000.biom")





#
#####################################  ALPHA AND BETA DIVERSITY  ###########################

############# Creating stacked area and stacked bars #################
#   ##### Define colors #####
nb.cols <- 40 #We know there are 15 ASVs here, so we can just keep them as such. 
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
mycolors #Pull color HEX codes from here. 
famcolsall = c( 'Acutalibacteraceae' = '#696969',
                "Akkermansiaceae" = "#A6CEE3",
                'Anaeroplasmataceae' = '#7FB5D5',
                'Anaerotignaceae' = '#6EAACF',
                'Anaerovoracaceae' = '#599DC8',
                "Bacteroidaceae" = '#5FA0CA',
                "Bifidobacteriaceae" = "#3386AE",
                'Borkfalkiaceae' = '#5BA2A2',
                'Burkholderiaceae' = '#72B29C',
                'CAG-272' = '#A5D981',
                'CAG-314' = '#A5D981',
                "CAG-508" = '#A5D981',
                'CAG-917' = '#A5D981',
                'Christensenellaceae' = '#91CE71',
                'Clostridiaceae_222000' = '#5EB54C',
                'Coprobacillaceae' = '#49AB3C',
                'Desulfovibrionaceae' = '#B15928',
                'Eggerthellaceae' = '#4F9F3B',
                'Enterobacteriaceae' = '#859D59',
                'Enterobacteriaceae_A' = '#859D59',
                'Enterococcaceae' = '#C1AF99',
                'Erysipelotrichaceae' = '#B89B74',
                'Gastranaerophilaceae' = '#DDB667', 
                'Lachnospiraceae' = '#D19B82', 
                'Lactobacillaceae' = '#F68383', 
                'Muribaculaceae' = '#A889C1',
                'Oscillospiraceae' = '#F48B55',
                "Oscillospiraceae_88309" = '#F48B55',
                'Nanosyncoccaceae' = '#FDBA67',
                'Peptostreptococcaceae' = '#F4892C',
                'Peptostreptococcaceae_256921' = '#F4892C',
                'Rikenellaceae' = '#CBB0CE',
                'Ruminococcaceae' = '#845DAA',
                'Staphylococcaceae' = '#977899',
                'Tannerellaceae' = '#B19A99',
                'Turicibacteraceae' = '#EAE499', 
                'Other' = '#000000')

#   ##### Stacked bars with microbiome package #####

devtools::install_github("microsud/microbiomeutilities")

#Tutorial from: https://microbiome.github.io/tutorials/Composition.html
# using microbiome package.
library(microbiome)
library(dplyr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(RColorBrewer)
library(microbiomeutilities) #version 1.0.17

#Get count of taxa
nrow(table(phyloseq::tax_table(pslach.2)[, "Genus"]))

#Decide on the color palette you want to call
display.brewer.all(colorblindFriendly=TRUE)
brewer.pal(n=8, "Set2")

# Make sure we use functions from correct package
transform <- microbiome::transform

# Merge rare taxa to speed up examples
pseq <- transform(psblanks, "compositional")
pseq1 <- aggregate_rare(pseq, level = "Family", detection = 0.001, prevalence = 0.001)

#Or if you want to show read counts instead of relative abundance:
psfam <- aggregate_rare(psblanks, level = "Family", detection = 0.01, prevalence = 0.01)

#pseq1 <- top_taxa(pseq, n=12)
rownames(table(phyloseq::tax_table(psfam)[, "Family"]))
nrow(table(phyloseq::tax_table(psfam)[, "Family"])) #Use this to see where your "Other" group falls. 

#Make sure grouping variable is recognized as a factor if you're faceting. 
pseq1@sam_data[["microb"]] <- factor(pseq1@sam_data[["microb"]], levels=c("LO","HI"))

#plot it. 
p <- plot_composition(psfam, taxonomic.level="Family", 
                      otu.sort= "abundance",
                      average_by = "sample_group") +
  scale_fill_manual(values=famcolsall) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=1),
        plot.title=element_text(size=16, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=14, face="bold"), 
        legend.text=element_text(size=13)) + 
  labs(x = "Blank Type", y = "Read count", title="psblanks averaged Stacked bar",  #Update these for other plots. 
       fill= "Family", subtitle="1U4U") + ggpubr::rotate_x_text(40) 
p
#print(p) +coord_flip()

ggsave("Stacked_bars/psblanksStackedFamily_abundance.pdf", p, units="in", height = 12, width=6)


#   ##### Making stacked area plots with microbiomeutilities #####
#      ##### By individual mouse. #####
ps.rel <- microbiome::transform(psexp3000hitreat, transform="compositional")
ps.rel@sam_data[["day_num"]] <- as.numeric(gsub("D", "", ps.rel@sam_data[["day"]]))


nrow(table(phyloseq::tax_table(ps.rel)[, "Family"]))
#Define colors if you need more than 12 total. 
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

plot_area(ps.rel, xvar="day_num", 
          level = "Family", 
          facet.by="mouse", 
          abund.thres = 0.01,
          prev.thres=0.1,
          fill.colors=mycolors, 
          nrow=4, ncol=5) + 
  ylab("Relative abundance") + 
  ggtitle("Stacked area of psexp3000hitreat") +
  scale_y_continuous(labels = scales::percent)




#

#      ##### Averaged by a grouping variable. #####
#         ##### Creating objects for hi and lo, all taxa #####
psexp3000hi.avg <- merge_samples(psexp3000hictrl, group="day", fun=mean)
psexp3000hi.avg@sam_data[["microb"]]="HI"
psexp3000hi.avg@sam_data[["day_num"]] <- as.numeric(gsub("D", "", rownames(sample_data(psexp3000hi.avg))))
sample_names(psexp3000hi.avg) <- paste0("HI", psexp3000hi.avg@sam_data[["day_num"]])

#View(sample_data(psexp3000hi.avg))

psexp3000lo.avg <- merge_samples(psexp3000loctrl, group="day", fun=mean)
psexp3000lo.avg@sam_data[["microb"]]="LO"
psexp3000lo.avg@sam_data[["day_num"]] <- as.numeric(gsub("D", "", rownames(sample_data(psexp3000lo.avg))))
sample_names(psexp3000lo.avg) <- paste0("LO", psexp3000lo.avg@sam_data[["day_num"]])
#Take the tree out to be able to merge the objects so that I can make a shared plot
# faceted by microbiome type. 
ps1_no_tree <- phyloseq(otu_table(psexp3000hi.avg), tax_table(psexp3000hi.avg), sample_data(psexp3000hi.avg))
ps2_no_tree <- phyloseq(otu_table(psexp3000lo.avg), tax_table(psexp3000lo.avg), sample_data(psexp3000lo.avg))
#Now merge the averaged object. 
pstreat.avg <- merge_phyloseq(ps1_no_tree, ps2_no_tree)
pstreat.avg <- prune_taxa(taxa_sums(pstreat.avg)>1, pstreat.avg)

#Make this into relative abundance. 
ps.rel.avg <- microbiome::transform(pstreat.avg, transform="compositional")

#Make sure your rows add to 1...
rowSums(otu_table(ps.rel.avg))

ps.rel.avg@sam_data[["microb"]] <- factor(ps.rel.avg@sam_data[["microb"]], levels=c("LO","HI"))
#

#         ##### Create your plot #####
#How many taxa total are there in the level of interest?
nrow(table(phyloseq::tax_table(ps.rel.avg)[, "Family"]))
#It'll combine the low abundance ones, but this is the max. 

#Define colors if you need more than 12 total. 
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

head(table(tax_table(ps.rel.avg)[, "Family"]))


#Note that x-variable has to be numeric, and facet.by is a required command for this function. 
p <- plot_area(ps.rel.avg, xvar="day_num", 
               level = "Family", 
               facet.by="microb", 
               abund.thres = 0.01,
               prev.thres=0.05,
               #otu.sort="abundance", #Can't figure out how to sort these. 
               fill.colors=famcolsall, 
               nrow=2, ncol=1) + 
  guides(col = guide_legend(ncol = 1, nrow = 23)) +
  xlab("Days post-infection") +
  ylab("Relative abundance") + 
  ggtitle("Stacked area all worm-treated animals") +
  scale_y_continuous(labels = scales::percent) + 
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=24, face="bold"), 
        axis.text=element_text(size=24), 
        legend.title=element_text(size=18, face="bold"), 
        legend.text=element_text(size=18)) +
  guides(fill=guide_legend(title="Family")) 

p

#
ggsave("Stacked_bars/psexp3000ctrl_StackedArea_MatchingColors.pdf", units = "in", height=14, width = 12)
#




############# Richness and evenness ################

#Summarize richness and evenness data
summary(dfprelo$evenness)
sd(dfprelo$evenness)

#Example of linear model analyzing alpha diversity at single timepoint
mod <- lm(Observed ~ microb + treatment + sex, data = dfpre)
AIC(mod)
summary(mod)


#Model richness and evenness with linear mixed effects models
lme.2 <- lme(PD ~ treatment + microb + day_num + treatment*day_num + microb*day_num, random=~1|mouse, data=dfexp3000)
anova(lme.2)



#







#

############# Beta diversity ################
#Load in packages. 
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(vegan)
library(ggbeeswarm)
library(plyr)
library(reshape2)
library(ggpmisc)
library(nlme)
library(lmerTest)
setwd("~/Documents/1U4U_16S/GG2_reclassification")

##### Compute distances between EACH TIMEPOINT AND PRETREATMENT - RESILIENCE #################

ps.rel <- microbiome::transform(psexp3000, transform="compositional")


#Make sure your objects are wht you want to be working with. 
pspre.rel
ps1.rel
ps10.rel
ps14.rel
ps28.rel
ps31.rel
ps41.rel

#Merge the phyloseq objects to get the preterat vs. each timepoint. 
pspre.1.rel <- merge_phyloseq(pspre.rel, ps1.rel)
pspre.10.rel <- merge_phyloseq(pspre.rel, ps10.rel)
pspre.14.rel <- merge_phyloseq(pspre.rel, ps14.rel)
pspre.28.rel <- merge_phyloseq(pspre.rel, ps28.rel)
pspre.31.rel <- merge_phyloseq(pspre.rel, ps31.rel)
pspre.41.rel <- merge_phyloseq(pspre.rel, ps41.rel)


#Get the beta diversity distances for each pairwise comparison back to 0. 
bc0.1 <- plasticity(pspre.1.rel, dist.method="bray", participant.col="mouse")
bc0.2 <- plasticity(pspre.10.rel, dist.method="bray", participant.col="mouse")
bc0.3 <- plasticity(pspre.14.rel, dist.method="bray", participant.col="mouse")
bc0.4 <- plasticity(pspre.28.rel, dist.method="bray", participant.col="mouse")
bc0.5 <- plasticity(pspre.31.rel, dist.method="bray", participant.col="mouse")
bc0.6 <- plasticity(pspre.41.rel, dist.method="bray", participant.col="mouse")

#   ##### Merging pairwise distance objects, all timepoints #####
bray.dm.time <- rbind(bc0.1, bc0.2, bc0.3, bc0.4, bc0.5, bc0.6) 

bray.dm.time$time_change <- (abs((bray.dm.time$day_num_2)-(bray.dm.time$day_num_1)))

#   ##### Now extract the GUnifrac distances. #####
library(GUniFrac)
ps.rel <- microbiome::transform(psexp3000, transform="compositional")
obj <- psexp3000

tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.

gunif <- GUniFrac(otu, tree, alpha=c(0,0.5,1))$unifracs #Generate the comparison matrix.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 

#Generate metadata that you can pull from down the line. 
obj <- psexp3000
dftmp <- data.frame(estimate_richness(obj, measures = "Observed"),
                    sum = sample_sums(obj), 
                    estimate_richness(obj, measures="Shannon"),
                    btools::estimate_pd(obj),    #This adds Faiths phylo div and spec richness
                    "row"=rownames(sample_data(obj)),
                    evenness=(estimate_richness(obj, measures="Shannon")/log(estimate_richness(obj, measures = "Observed"))),
                    "cage_day"=paste0(sample_data(obj)$cage,sample_data(obj)$day))  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
exp_dat<- cbind(sample_data(psexp3000),dftmp)
head(exp_dat)

exp_dat$microb <- factor(exp_dat$microb, levels=c("LO", "HI"))

#Pull the distances. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#row=rownames(sample_data(psexp3000))
sd = exp_dat %>%    #Need to change this when doing a different physeq object
  select(row, div_day, microb, treatment, day, mouse) %>%
  mutate_if(is.factor,as.character)
colnames(sd) = c("Var1", "microb_day", "microb", "treatment", "day","mouse")
GU50.sd = left_join(GU50.m, sd, by = "Var1")
colnames(sd) = c("Var2", "microb_day", "microb", "treatment", "day","mouse")
GU50.sd = left_join(GU50.sd, sd, by = "Var2")
GU50.sd$MT <- paste(GU50.sd$microb.x, GU50.sd$treatment.x, sep="_") #add another useful column for plotting. 

#Pull the comparisons of the individual change over all timepoints. 
GU50.sd.1 <- GU50.sd[GU50.sd$mouse.x==GU50.sd$mouse.y & 
                       GU50.sd$day.x!=GU50.sd$day.y &
                       GU50.sd$day.x=="D-5" ,]
GU50.sd.1 <- GU50.sd.1 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1)
GU50.sd.2 <- GU50.sd[GU50.sd$mouse.x==GU50.sd$mouse.y & 
                       GU50.sd$day.x!=GU50.sd$day.y &
                       GU50.sd$day.y=="D-5" ,]
GU50.sd.final <- rbind(GU50.sd.1, GU50.sd.2)

#Subset to merge with other distance data. 
GU50.sub <- GU50.sd.final %>%
  select(Var1, Var2, value) 
GU50.sub <- GU50.sub %>%
  dplyr::rename(S1=Var1, 
                S2=Var2, 
                gunif=value)


#Put it all together and then add more metadata columns. 
dist.dm.time <- bray.dm.time
#ok now merge the GUnifrac distances with the others. 
dist.dm.time <- left_join(dist.dm.time, GU50.sub, by="S1") #Pull in the GUnifrac values. 
colnames(dist.dm.time) <- make.unique(names(dist.dm.time))
dist.dm.time$MT <- paste(dist.dm.time$microb_1, dist.dm.time$treatment_1, sep="_")
dist.dm.time$MT <- factor(dist.dm.time$MT, levels=c("LO_control","LO_worm", "HI_control", "HI_worm"))
dist.dm.time$microb_1 <- factor(dist.dm.time$microb_1, levels=c("LO", "HI"))
dist.dm.time$div_day_1 <- factor(dist.dm.time$div_day_1, levels=c('LO_D-5', 'LO_D1', 'LO_D10', 'LO_D14', 'LO_D28', 
                                                                  'LO_D31', 'LO_D41', 'HI_D-5', 'HI_D1', 'HI_D10', 
                                                                  'HI_D14', 'HI_D28', 'HI_D31', 'HI_D41'))
dist.dm.time$MDT <- paste(dist.dm.time$microb_1, dist.dm.time$day_1, dist.dm.time$treatment_1, sep="_")
dist.dm.time$MDT <- factor(dist.dm.time$MDT, 
                           levels=c('LO_D-5_worm', 'LO_D1_worm', 'LO_D10_worm', 'LO_D14_worm', 'LO_D28_worm','LO_D31_worm', 'LO_D41_worm',
                                    'HI_D-5_worm', 'HI_D1_worm', 'HI_D10_worm', 'HI_D14_worm', 'HI_D28_worm', 'HI_D31_worm', 'HI_D41_worm', 
                                    'LO_D-5_control', 'LO_D1_control', 'LO_D10_control', 'LO_D14_control', 'LO_D28_control', 'LO_D31_control', 'LO_D41_control', 
                                    'HI_D-5_control', 'HI_D1_control', 'HI_D10_control', 'HI_D14_control', 'HI_D28_control', 'HI_D31_control', 'HI_D41_control'))
dist.dm.time$treatment_1 <- factor(dist.dm.time$treatment_1, levels=c("worm", "control"))
dist.dm.time$day_num_1 <- as.numeric(dist.dm.time$day_num_1)
dist.dm.time$div_day_1 <- factor(dist.dm.time$div_day_1, levels=c("LO_D1",  "LO_D10", "LO_D14", "LO_D28", "LO_D31", "LO_D41", 
                                                                  "HI_D1", "HI_D10", "HI_D14", "HI_D28", "HI_D31", "HI_D41"))


#  


#   ##### Subsetting whole df into various variables #####
dist.dm.time.treat <- dist.dm.time[dist.dm.time$treatment_1 =="worm",]
dist.dm.time.ctrl <- dist.dm.time[dist.dm.time$treatment_1 =="control",]

dist.dm.time.hi <- dist.dm.time[dist.dm.time$microb_1 =="HI",]
dist.dm.time.lo <- dist.dm.time[dist.dm.time$microb_1 =="LO",]

#All timepoints LO and HI, by microbiome diversity and treatment
lotreat <- dist.dm.time.treat[dist.dm.time.treat$microb_1=="LO",]
hitreat <- dist.dm.time.treat[dist.dm.time.treat$microb_1=="HI",]
losham <- dist.dm.time.ctrl[dist.dm.time.ctrl$microb_1=="LO",]
hisham <- dist.dm.time.ctrl[dist.dm.time.ctrl$microb_1=="HI",]

#Note that these are only worm-treated. 
#Subset just to pretreatment, infected, and cleared timepoints for worm-treated 
dist.dm.imp.treat <- subset(dist.dm.time.treat, day_1 %in% c("D14", "D41"))
dist.dm.pre14.treat <- dist.dm.imp.treat[dist.dm.imp.treat$day_1 == "D14",]
dist.dm.pre41.treat <- dist.dm.imp.treat[dist.dm.imp.treat$day_1 == "D41",]
#
#Do the same for sham-treated
dist.dm.imp.ctrl <- subset(dist.dm.time.ctrl, day_1 %in% c("D14", "D41"))
dist.dm.pre14.ctrl <- dist.dm.imp.ctrl[dist.dm.imp.ctrl$day_1 == "D14",]
dist.dm.pre41.ctrl <- dist.dm.imp.ctrl[dist.dm.imp.ctrl$day_1 == "D41",]

#Now only LOW animals.
dist.dm.imp.lo <- subset(dist.dm.time, day_1 %in% c("D14", "D41") & microb_1 %in% c("LO"))
dist.dm.pre14.lo <- dist.dm.imp.lo[dist.dm.imp.lo$day_1 == "D14",]
dist.dm.pre41.lo <- dist.dm.imp.lo[dist.dm.imp.lo$day_1 == "D41",]

#Now only HI animals
dist.dm.imp.hi <- subset(dist.dm.time, day_1 %in% c("D14", "D41") & microb_1 %in% c("HI"))
dist.dm.pre14.hi <- dist.dm.imp.hi[dist.dm.imp.hi$day_1 == "D14",]
dist.dm.pre41.hi <- dist.dm.imp.hi[dist.dm.imp.hi$day_1 == "D41",]

#


#   ##### Plotting boxplots of pre to infection, and pre to clearance. #####
dist.dm.imp.hi
dist.dm.imp.lo
dist.dm.imp.treat
dist.dm.imp.ctrl
dist.dm.imp.treat$microb_day <- paste0(dist.dm.imp.treat$microb_1, dist.dm.imp.treat$day_1)
dist.dm.imp.ctrl$microb_day <- paste0(dist.dm.imp.ctrl$microb_1, dist.dm.imp.ctrl$day_1)


loimpcols <- c('#FDF8FF' , '#A11FE5', '#300049')
hiimpcols <- c("#FFFCF3" , "#DCAE1C" , "#644C00")

loimpellipse <- c('#c2b6cf' , '#A11FE5', '#300049')
hiimpellipse <- c("#c9c4b3" , "#DCAE1C" , "#644C00")


#Comparing worm-treated between pretreat to infection, then pretreat to clearance. 
ggplot(dist.dm.imp.treat, aes(x=microb_1, y=bray)) +
  geom_boxplot(aes(colour=div_day_1), outlier.shape=NA) + 
  geom_beeswarm(aes(shape=treatment_1, fill=div_day_1), cex=4, size=8, alpha=0.9) +
  #geom_point(aes(shape=sex, fill=microb), position=position_jitter(0.1), size=10, alpha=0.9) +
  scale_colour_manual(values=c( "#bf9001", "#644C00", "#7e5ba0",'#300049' )) +
  scale_fill_manual(values=c( "#bf9001", "#644C00", "#7e5ba0",'#300049')) +
  scale_shape_manual(values=c(21)) +
  ylim(min(dist.dm.imp.treat$bray), max(dist.dm.imp.treat$bray)+1*sd(dist.dm.imp.treat$bray)) +
  ggtitle("Pre to inf and clearance Δ bray, worm-treated") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position="none") +
  facet_wrap(~day_1)


#All hi, comparing worm- vs. sham-treated and facted by day
ggplot(dist.dm.imp.hi, aes(x=treatment_1, y=bray)) +
  geom_boxplot(aes(colour=day_1), outlier.shape=NA) + 
  geom_beeswarm(aes(shape=treatment_1, fill=day_1), cex=4, size=8, alpha=0.9) +
  scale_colour_manual(values=c("#bf9001", "#644C00")) +
  scale_fill_manual(values=c("#bf9001", "#644C00")) +
  scale_shape_manual(values=c(24,21)) +
  ylim(min(dist.dm.imp.hi$bray), max(dist.dm.imp.hi$bray)+1*sd(dist.dm.imp.hi$bray)) +
  ggtitle("Pre to inf and clearance Δ bray, HI") +
  theme_bw() +
  theme(text=element_text(size=24), legend.position="none") +
  facet_wrap(~day_1)


#
#   ##### Plotting line plots of change since start #####

mem.hitreat <- lmer(bray~ day_num_1 + (1|mouse_1), data=dist.dm.time.hi)
mem.lotreat <- lmer(bray~ day_num_1 + (1|mouse_1), data=dist.dm.time.lo)
mem.losham <- lmer(bray~ day_num_1 + (1|mouse_1), data=losham)
mem.hisham <- lmer(bray~ day_num_1 + (1|mouse_1), data=hisham)

fixef(mem.lotreat)
fixef(mem.hitreat)
fixef(mem.losham)
fixef(mem.hisham)


p12 = ggplot(dist.dm.time, aes(x=day_num_1, y=bray)) + 
  geom_beeswarm(aes(fill=MDT, shape=microb_1), dodge.width=3, alpha=0.8, size=5, na.rm = TRUE, cex=1.1) +  
  ylab("Individual dissimilarity score") + 
  xlab("Days post-infection") + 
  scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=allcols1) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text.y=element_text(size=20)) +
  #stat_poly_eq(use_label(c("eq", "R2"))) +  
  facet_wrap(~microb_1) +
  ggtitle("Bray-Curtis Change since pretreat")
#ok now make sure you have your lmer results. Use "Estimate" column.
p12allcols <- p12 + 
  geom_abline(intercept=data.frame(fixef(mem.lotreat))[1,1], slope=data.frame(fixef(mem.lotreat))[2,1], colour="#A11FE5") + #Outcome of lotreat only model
  geom_abline(intercept=data.frame(fixef(mem.hitreat))[1,1], slope=data.frame(fixef(mem.hitreat))[2,1], colour="#DCAE1C") + #Outcome of hitreat
  geom_abline(intercept=data.frame(fixef(mem.losham))[1,1], slope=data.frame(fixef(mem.losham))[2,1], colour="#383838") + #Outcome of losham
  geom_abline(intercept=data.frame(fixef(mem.hisham))[1,1], slope=data.frame(fixef(mem.hisham))[2,1], colour="#292929")  #Outcome of hisham
p12allcols




################## Running PERMANOVAs on GUniFrac distances and Bray-Curtis #####

#   ##### Running GUniFrac PERMANOVA with adonis3 #####
metadata <- as(sample_data(psexp3000), "data.frame")  
perms <- with(metadata, how(nperm = 5000)) #Note increasing number of permutations makes results more consistent.
mod <-  adonis3(d50 ~ microb +lane + day_num + treatment, data=metadata)
mod$aov.tab #View your p-values. 


#   ##### Running regular vegan PERMANOVA with adonis2 for other metrics #####
metadata <- as(sample_data(psprehi), "data.frame")  
perms <- with(metadata, how(nperm = 5000)) #Note increasing number of permutations makes results more consistent.
mod <- adonis2(phyloseq::distance(psprehi, method="bray") ~  lane + treatment, data=metadata)
mod


#



################################  Code from other workpaces, with workspace loading and colors defined for each  #############################
############ Looking for differences in body mass over time ###############
setwd("~/Documents/1U4U_16S/Rerun_metadata/")

dat <- read.csv("1U4U_body_mass.csv")
dat$div_day <- paste(dat$microb, dat$day_num, sep="_")
dat$div_day <- factor(dat$div_day, levels=c("lo_-5" , "lo_14", "lo_41" , "hi_-5", "hi_14", "hi_41"))
dat$MTD <- paste(dat$div_day, dat$treatment, sep="_")
dat$MTD <- factor(dat$MTD, levels=c("lo_-5_worm","lo_14_worm", "lo_41_worm",
                                    "lo_-5_control", "lo_14_control", "lo_41_control" , 
                                    "hi_-5_worm", "hi_14_worm", "hi_41_worm",  
                                    "hi_-5_control","hi_14_control", "hi_41_control"))
dat$microb <- factor(dat$microb, levels=c("lo", "hi"))
dat$day <- paste0("D", dat$day_num)
colnames(dat)

divdaycols<- c('#FDF8FF' , '#A11FE5' , '#300049', 
               "#FFFCF3" , "#DCAE1C" , "#644C00")
MTDcols <- c('#FDF8FF' , '#A11FE5' , '#300049', 
             '#FFFFFF', '#FFFFFF', '#FFFFFF', 
             "#FFFCF3" , "#DCAE1C" , "#644C00", 
             '#FFFFFF', '#FFFFFF', '#FFFFFF')
MTDcols1 <- c('#ECC9FE' , '#A11FE5' , '#300049', 
              '#ECC9FE' , '#A11FE5' , '#300049',
              "#FFEFBA" , "#DCAE1C" , "#644C00", 
              "#FFEFBA" , "#DCAE1C" , "#644C00" )


ggplot(data=dat, aes(x=day_num, y = mass)) +
  geom_boxplot(aes(color=MTD), width=4, outlier.shape=NA) +
  geom_beeswarm(aes(fill=MTD, shape=microb), dodge.width=4, alpha=0.8, size=5, na.rm = TRUE, cex=1.1) +  
  scale_color_manual(values=MTDcols1) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=MTDcols) +
  theme_bw() +
  facet_wrap(~microb) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text.y=element_text(size=20)) +
  labs(x="Sampling Day", y="Body mass (g)", title="Body Mass All Animals 1U4U")


ggsave("~/Documents/1U4U_16S/GG2_reclassification/bodymassBoxplots_psimptimepoints.pdf", units="in", height=5, width=9)


#Modelling the outcomes that I'm looking at. 
datlo <- dat[dat$microb=="lo",]
dathi <- dat[dat$microb=="hi",]

mem.bod <- lme(mass ~ treatment + day_num, random=~1|mouse, data=datlo)
mem.bod
anova(mem.bod)

emmeans(mem.bod, pairwise ~ day_num)

############ Looking at microbial load (copies of 16S) from qPCR ###############

setwd("~/Documents/1U4U_qPCR")

#Read in the whole df. 

allqpcr <- read.csv("./1U4U_Full16S_metadata_ByOrigID_wBactQuant_wAll18684.csv") #original data from Zac. 
newqpcr <- read.csv("map_1U4U_Full16S_metadata_ByOrigID_wBactQuant.csv") #Pull the read counts. 
dfps1 <- read.csv("/Users/mdoolin/Documents/1U4U_16S/Rerun_metadata/ps1_metadata.csv")
div_day <- paste(dfps1$microb, sep="_", dfps1$day)
reads <- cbind("GnomexID"=dfps1$GNomexID, "ASVs"=dfps1$Observed,
               "count"=dfps1$sum, div_day)

#Subset for only samples with more than 3000 reads. This is what's in the rarefied data.
qpcr3000 <- allqpcr.1[allqpcr.1$count>3000,]

qpcr3000$div_day <- factor(qpcr3000$div_day, levels=c("LO_D-5", "LO_D1",  "LO_D10","LO_D14" ,"LO_D28", "LO_D31", "LO_D41", 
                                                      "HI_D-5", "HI_D1", "HI_D10","HI_D14" ,"HI_D28", "HI_D31" ,"HI_D41"))
qpcr3000$day_num <-as.numeric(gsub("D", "", qpcr3000$day))
qpcr3000$microb <- factor(qpcr3000$microb, levels=c("LO", "HI"))

qpcr.imp <- subset(qpcr3000, day %in% c("D-5", "D14", "D41"))
qpcr.imp$MTD <- paste(qpcr.imp$microb, qpcr.imp$day_num, qpcr.imp$treatment, sep="_")
qpcr.imp$MTD <- factor(qpcr.imp$MTD, levels=c("LO_-5_worm","LO_14_worm", "LO_41_worm",
                                              "LO_-5_control", "LO_14_control", "LO_41_control" , 
                                              "HI_-5_worm", "HI_14_worm", "HI_41_worm",  
                                              "HI_-5_control","HI_14_control", "HI_41_control"))

MTDcols <- c('#FDF8FF' , '#A11FE5' , '#300049', #For use when defining "fill" with 3 timepoints of treatment and controls. 
             '#FFFFFF', '#FFFFFF', '#FFFFFF', 
             "#FFFCF3" , "#DCAE1C" , "#644C00", 
             '#FFFFFF', '#FFFFFF', '#FFFFFF')
MTDcols1 <- c('#ECC9FE' , '#A11FE5' , '#300049', #For use when defining colors with 3 timepoints of treatment and controls. 
              '#ECC9FE' , '#A11FE5' , '#300049',
              "#FFEFBA" , "#DCAE1C" , "#644C00", 
              "#FFEFBA" , "#DCAE1C" , "#644C00" )



impplot <- ggplot(qpcr.imp, aes(x=day_num, y=copies_16S)) +
  geom_boxplot(aes(color=MTD), width=4, outlier.shape=NA) +
  geom_beeswarm(aes(fill=MTD, shape=microb), dodge.width=4, alpha=0.8, size=5, na.rm = TRUE, cex=1.1) +  
  scale_color_manual(values=MTDcols1) +
  scale_shape_manual(values=c(21,22)) +
  scale_fill_manual(values=MTDcols) +
  theme_bw() +
  facet_wrap(~microb) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text.y=element_text(size=20)) +
  labs(x="Sampling Day", y="Number 16S copies", title="16S copies by qPCR, psimp")
impplot

ggsave("psimpBoxplots_Copies16S.pdf", units="in", height=5, width=9)



# Running the models on the various days. 
library(nlme)
mem.load <- lme(copies_16S ~ microb + treatment + day_num + treatment * day_num+ microb*day_num, random = ~1|mouse, data=qpcr.imp)
anova(mem.load)

##            numDF denDF  F-value p-value
##(Intercept)     1    90 334.6539  <.0001
##microb          1    57   5.0828  0.0280
##treatment       1    57   0.0004  0.9841
##day_num         1    90   3.9500  0.0499

emmeans(mem.load, pairwise ~ microb)
#day is significant because D14 is lower than the other two days. 
#microb is significant because LO are lower than HI on D14 regardless of treatment. 

#Now just looking at the LO and HI individually to see if there is an effect of day or treatment. 
qpcr.implo <- qpcr.imp[qpcr.imp$microb=="LO",]
qpcr.imphi <- qpcr.imp[qpcr.imp$microb=="HI",]

mem.load <- lme(copies_16S ~ treatment + day_num, random = ~1|mouse, data=qpcr.imphi)
anova(mem.load)

##For LO animals only
##            numDF denDF   F-value p-value
##(Intercept)     1    63 193.19013  <.0001
##treatment       1    36   0.00383  0.9510
##day_num         1    63   1.60326  0.2101

##For HI animals only
##           numDF denDF   F-value p-value
##(Intercept)     1    26 139.20454  <.0001
##treatment       1    20   0.01262  0.9117
##day_num         1    26   2.82163  0.1050


#

################################################## 16S Data ############################################
############ Generating the lineplots for change in beta diversity plots and models for those data. ######

#Load in packages. 
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(vegan)
library(ggbeeswarm)
library(plyr)
library(reshape2)
library(ggpmisc)
setwd("~/Documents/1U4U_16S/GG2_reclassification")

psexp3000 <- readRDS("psexp3000wo18RLP.rds")
psexp3000@sam_data[["day_num"]] <- as.numeric(gsub("D", "", psexp3000@sam_data[["day"]]))
ps.rel <- microbiome::transform(psexp3000, transform="compositional")

rename=dplyr::rename
#


# Define color palettes. #
#Colors, in order
locols <- c('#d6a3f0' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
hicols <- c("#f5e1a6" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")

allcols <- c('#d6a3f0' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049',
             "#f5e1a6" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")

#Colors to use to define all control samples as white, with variable being "MDT". 
allcols1 <- c('#d6a3f0' , '#c275eb' , '#A11FE5' , '#7E05BE' , '#5C008D', '#300049',
              "#f5e1a6" , "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00", 
              "#FFFFFF", "#FFFFFF", "#FFFFFF",  "#FFFFFF", "#FFFFFF", "#FFFFFF",
              "#FFFFFF", "#FFFFFF", "#FFFFFF",  "#FFFFFF", "#FFFFFF", "#FFFFFF")


#    Make sure your objects are what you want to be working with. Subset to individual timepoints. # 
pspre.rel <- subset_samples(ps.rel, day_num==-5)
ps1.rel <- subset_samples(ps.rel, day_num==1)
ps10.rel <- subset_samples(ps.rel, day_num==10)
ps14.rel <- subset_samples(ps.rel, day_num==14)
ps28.rel <- subset_samples(ps.rel, day_num==28)
ps31.rel <- subset_samples(ps.rel, day_num==31)
ps41.rel <- subset_samples(ps.rel, day_num==41)

#Merge the phyloseq objects to get the preterat vs. each timepoint. 
pspre.1.rel <- merge_phyloseq(pspre.rel, ps1.rel)
pspre.10.rel <- merge_phyloseq(pspre.rel, ps10.rel)
pspre.14.rel <- merge_phyloseq(pspre.rel, ps14.rel)
pspre.28.rel <- merge_phyloseq(pspre.rel, ps28.rel)
pspre.31.rel <- merge_phyloseq(pspre.rel, ps31.rel)
pspre.41.rel <- merge_phyloseq(pspre.rel, ps41.rel)

#Get the beta diversity distances for each pairwise comparison back to 0. 
bc0.1 <- plasticity(pspre.1.rel, dist.method="bray", participant.col="mouse")
bc0.2 <- plasticity(pspre.10.rel, dist.method="bray", participant.col="mouse")
bc0.3 <- plasticity(pspre.14.rel, dist.method="bray", participant.col="mouse")
bc0.4 <- plasticity(pspre.28.rel, dist.method="bray", participant.col="mouse")
bc0.5 <- plasticity(pspre.31.rel, dist.method="bray", participant.col="mouse")
bc0.6 <- plasticity(pspre.41.rel, dist.method="bray", participant.col="mouse")

U0.1 <- plasticity(pspre.1.rel, dist.method="unifrac", participant.col="mouse")
U0.2 <- plasticity(pspre.10.rel, dist.method="unifrac", participant.col="mouse")
U0.3 <- plasticity(pspre.14.rel, dist.method="unifrac", participant.col="mouse")
U0.4 <- plasticity(pspre.28.rel, dist.method="unifrac", participant.col="mouse")
U0.5 <- plasticity(pspre.31.rel, dist.method="unifrac", participant.col="mouse")
U0.6 <- plasticity(pspre.41.rel, dist.method="unifrac", participant.col="mouse")

WU0.1 <- plasticity(pspre.1.rel, dist.method="wunifrac", participant.col="mouse")
WU0.2 <- plasticity(pspre.10.rel, dist.method="wunifrac", participant.col="mouse")
WU0.3 <- plasticity(pspre.14.rel, dist.method="wunifrac", participant.col="mouse")
WU0.4 <- plasticity(pspre.28.rel, dist.method="wunifrac", participant.col="mouse")
WU0.5 <- plasticity(pspre.31.rel, dist.method="wunifrac", participant.col="mouse")
WU0.6 <- plasticity(pspre.41.rel, dist.method="wunifrac", participant.col="mouse")

# Merging pairwise distance objects, all timepoints #
bray.dm.time <- rbind(bc0.1, bc0.2, bc0.3, bc0.4, bc0.5, bc0.6) 
jaccard.dm.time <- rbind(j0.1, j0.2, j0.3, j0.4, j0.5, j0.6)
unifrac.dm.time <- rbind(U0.1, U0.2, U0.3, U0.4, U0.5, U0.6)
wunifrac.dm.time <- rbind(WU0.1, WU0.2, WU0.3, WU0.4, WU0.5, WU0.6) 

bray.dm.time$time_change <- (abs((bray.dm.time$day_num_2)-(bray.dm.time$day_num_1)))
jaccard.dm.time$time_change <- (abs((jaccard.dm.time$day_num_2)-(jaccard.dm.time$day_num_1)))
unifrac.dm.time$time_change <- (abs((unifrac.dm.time$day_num_2)-(unifrac.dm.time$day_num_1)))
wunifrac.dm.time$time_change <- (abs((wunifrac.dm.time$day_num_2)-(wunifrac.dm.time$day_num_1)))


#      Ok now extract the GUnifrac distances. #
library(GUniFrac)
#ps.rel <- microbiome::transform(psexp3000, transform="compositional")
obj <- psexp3000

tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.

gunif <- GUniFrac(otu, tree, alpha=c(0,0.5,1))$unifracs #Generate the comparison matrix.
#d100 <- gunif[, , "d_1"]
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
#d75 <- gunif[, , "d_0.75"]

#     Generate metadata that you can pull from down the line. 
obj <- psexp3000
dftmp <- data.frame(estimate_richness(obj, measures = "Observed"),
                    sum = sample_sums(obj), 
                    estimate_richness(obj, measures="Shannon"),
                    btools::estimate_pd(obj),    #This adds Faiths phylo div and spec richness
                    "row"=rownames(sample_data(obj)),
                    evenness=(estimate_richness(obj, measures="Shannon")/log(estimate_richness(obj, measures = "Observed"))),
                    "cage_day"=paste0(sample_data(obj)$cage,sample_data(obj)$day))  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
exp_dat<- cbind(sample_data(psexp3000),dftmp)
head(exp_dat)

exp_dat$microb <- factor(exp_dat$microb, levels=c("LO", "HI"))

#We only really care about the 50% weighting to compare to unweighted and weighted unifrac. 
#Pull the distances. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#row=rownames(sample_data(psexp3000))
sd = exp_dat %>%    #Need to change this when doing a different physeq object
  select(row, div_day, microb, treatment, day, mouse) %>%
  mutate_if(is.factor,as.character)
colnames(sd) = c("Var1", "microb_day", "microb", "treatment", "day","mouse")
GU50.sd = left_join(GU50.m, sd, by = "Var1")
colnames(sd) = c("Var2", "microb_day", "microb", "treatment", "day","mouse")
GU50.sd = left_join(GU50.sd, sd, by = "Var2")
GU50.sd$MT <- paste(GU50.sd$microb.x, GU50.sd$treatment.x, sep="_") #add another useful column for plotting. 

#Pull the comparisons of the individual change over all timepoints. 
GU50.sd.1 <- GU50.sd[GU50.sd$mouse.x==GU50.sd$mouse.y & 
                       GU50.sd$day.x!=GU50.sd$day.y &
                       GU50.sd$day.x=="D-5" ,]
GU50.sd.1 <- GU50.sd.1 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1)
GU50.sd.2 <- GU50.sd[GU50.sd$mouse.x==GU50.sd$mouse.y & 
                       GU50.sd$day.x!=GU50.sd$day.y &
                       GU50.sd$day.y=="D-5" ,]
GU50.sd.final <- rbind(GU50.sd.1, GU50.sd.2)

#Subset to merge with other distance data. 
GU50.sub <- GU50.sd.final %>%
  select(Var1, Var2, value) 
GU50.sub <- GU50.sub %>%
  dplyr::rename(S1=Var1, 
                S2=Var2, 
                gunif=value)


#     Put it all together and then add more metadata columns. 
dist.dm.time <- cbind(bray.dm.time, jaccard.dm.time$jaccard, unifrac.dm.time$unifrac, wunifrac.dm.time$wunifrac)
colnames(dist.dm.time)[which(names(dist.dm.time)=="jaccard.dm.time$jaccard")] <- "jaccard"
colnames(dist.dm.time)[which(names(dist.dm.time)=="unifrac.dm.time$unifrac")] <- "unifrac"
colnames(dist.dm.time)[which(names(dist.dm.time)=="wunifrac.dm.time$wunifrac")] <- "wunifrac"

#ok now merge the GUnifrac distances with the others. 
dist.dm.time <- left_join(dist.dm.time, GU50.sub, by="S1") #Pull in the GUnifrac values. 
colnames(dist.dm.time) <- make.unique(names(dist.dm.time))
dist.dm.time$MT <- paste(dist.dm.time$microb_1, dist.dm.time$treatment_1, sep="_")
dist.dm.time$MT <- factor(dist.dm.time$MT, levels=c("LO_control","LO_worm", "HI_control", "HI_worm"))
dist.dm.time$microb_1 <- factor(dist.dm.time$microb_1, levels=c("LO", "HI"))
dist.dm.time$div_day_1 <- factor(dist.dm.time$div_day_1, levels=c('LO_D-5', 'LO_D1', 'LO_D10', 'LO_D14', 'LO_D28', 
                                                                  'LO_D31', 'LO_D41', 'HI_D-5', 'HI_D1', 'HI_D10', 
                                                                  'HI_D14', 'HI_D28', 'HI_D31', 'HI_D41'))
dist.dm.time$MDT <- paste(dist.dm.time$microb_1, dist.dm.time$day_1, dist.dm.time$treatment_1, sep="_")
dist.dm.time$MDT <- factor(dist.dm.time$MDT, 
                           levels=c('LO_D-5_worm', 'LO_D1_worm', 'LO_D10_worm', 'LO_D14_worm', 'LO_D28_worm','LO_D31_worm', 'LO_D41_worm',
                                    'HI_D-5_worm', 'HI_D1_worm', 'HI_D10_worm', 'HI_D14_worm', 'HI_D28_worm', 'HI_D31_worm', 'HI_D41_worm', 
                                    'LO_D-5_control', 'LO_D1_control', 'LO_D10_control', 'LO_D14_control', 'LO_D28_control', 'LO_D31_control', 'LO_D41_control', 
                                    'HI_D-5_control', 'HI_D1_control', 'HI_D10_control', 'HI_D14_control', 'HI_D28_control', 'HI_D31_control', 'HI_D41_control'))
dist.dm.time$treatment_1 <- factor(dist.dm.time$treatment_1, levels=c("worm", "control"))
dist.dm.time$day_num_1 <- as.numeric(dist.dm.time$day_num_1)
dist.dm.time$div_day_1 <- factor(dist.dm.time$div_day_1, levels=c("LO_D1",  "LO_D10", "LO_D14", "LO_D28", "LO_D31", "LO_D41", 
                                                                  "HI_D1", "HI_D10", "HI_D14", "HI_D28", "HI_D31", "HI_D41"))
#  


#      Subsetting whole df into various variables #
dist.dm.time.treat <- dist.dm.time[dist.dm.time$treatment_1 =="worm",]
dist.dm.time.ctrl <- dist.dm.time[dist.dm.time$treatment_1 =="control",]

dist.dm.time.hi <- dist.dm.time[dist.dm.time$microb_1 =="HI",]
dist.dm.time.lo <- dist.dm.time[dist.dm.time$microb_1 =="LO",]

#All timepoints LO and HI, by microbiome diversity and treatment
lotreat <- dist.dm.time.treat[dist.dm.time.treat$microb_1=="LO",]
hitreat <- dist.dm.time.treat[dist.dm.time.treat$microb_1=="HI",]
losham <- dist.dm.time.ctrl[dist.dm.time.ctrl$microb_1=="LO",]
hisham <- dist.dm.time.ctrl[dist.dm.time.ctrl$microb_1=="HI",]


#     IMPORTANT MODELS and PLOTS: Bray-Curtis #
#Generate your linear mixed effects model outcomes including individual as a random effect
# and time (numeric) as the fixed effect. Record the model outcomes elsewhere. 
mem.hitreat <- lmer(bray~ day_num_1 + (1|mouse_1), data=dist.dm.time.hi)
mem.lotreat <- lmer(bray~ day_num_1 + (1|mouse_1), data=dist.dm.time.lo)
mem.losham <- lmer(bray~ day_num_1 + (1|mouse_1), data=losham)
mem.hisham <- lmer(bray~ day_num_1 + (1|mouse_1), data=hisham)
#You can look at the fixed effects here. They're pulled again to generate the correct lines for plots
# within the plotting code. 
fixef(mem.lotreat)
fixef(mem.hitreat)
fixef(mem.losham)
fixef(mem.hisham)

#Plot the data. 
p12 = ggplot(dist.dm.time, aes(x=day_num_1, y=bray)) + 
  geom_beeswarm(aes(fill=MDT, shape=microb_1), dodge.width=3, alpha=0.8, size=5, na.rm = TRUE, cex=1.1) +  
  #geom_point(aes(shape=treatment_1, fill=MDT), size=2, alpha=0.8, position=position_dodge(width=4)) +
  #theme_biome_utils() + 
  ylab("Individual dissimilarity score") + 
  xlab("Days post-infection") + 
  scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=allcols1) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text.y=element_text(size=20)) +
  facet_wrap(~microb_1) +
  ggtitle("Bray-Curtis Change since pretreat")
#ok now make sure you have your lmer results. Use "Estimate" column.
p12allcols <- p12 + 
  geom_abline(intercept=data.frame(fixef(mem.lotreat))[1,1], slope=data.frame(fixef(mem.lotreat))[2,1], colour="#A11FE5") + #Outcome of lotreat only model
  geom_abline(intercept=data.frame(fixef(mem.hitreat))[1,1], slope=data.frame(fixef(mem.hitreat))[2,1], colour="#DCAE1C") + #Outcome of hitreat
  geom_abline(intercept=data.frame(fixef(mem.losham))[1,1], slope=data.frame(fixef(mem.losham))[2,1], colour="#6C696E") + #Outcome of losham
  geom_abline(intercept=data.frame(fixef(mem.hisham))[1,1], slope=data.frame(fixef(mem.hisham))[2,1], colour="#6C696e")   #Outcome of hisham
p12allcols


#Then do the same for the others. 

# Arrange the plots with all points in one #

g <- arrangeGrob(p12allcols, p13allcols, p14allcols, p15allcols, ncol=2) #This is so that you can save it manually. 
grid.arrange(p12allcols, p13allcols, p14allcols, p15allcols, ncol=2) #This is so you can look at it. 

ggsave("BetaDiv/Lineplots/psexp3000_LinePlotsBrayWUUnifracGUnif_allcols.pdf", g, units="in", width=28, height=14)

#


#      Modeling for the whole dataset and partial datasets. #
#For the plots, I generated simple lines just for the groups within the whole dataset for which I wanted to show
# trends. 
hist(dist.dm.time$bray)
mem.bray <- lme(bray ~ microb_1*day_num_1 + treatment_1*day_num_1 + 
                  microb_1*treatment_1, random=~1|mouse_1, data = dist.dm.time)
summary(mem.bray)


#Only looking at the worm-treated animals:
mem.1 <- lme(gunif ~ day_num_1 * microb_1, random=~1|mouse_1, data=dist.dm.time.treat)
summary(mem.1)



#Can show pairwise post-hoc within each effect using emmeans. 
library(emmeans)
emmeans(mod, list(pairwise ~ day_num_1:microb_1),)

#






############ Alpha and beta diversity plotting ######
#   ##### Richness/Evenness Scatterplot #####
setwd("~/Documents/1U4U_16S/GG2_reclassification/")

#load libraries
library(easypackages)
libraries("phyloseq", "ggplot2", "dplyr", "ggpubr", "gridExtra", "microbiome", 
          "scales", "microbiomeutilities", "vegan", "GUniFrac")
#vegan version 2.6.4

#Load in main phyloseq objects. 
psexp3000 <- readRDS("R_objects/psexp3000_gg2.rds")



obj <- psexp3000

dftmp <- data.frame(estimate_richness(obj, measures = "Observed"),
                    sum = sample_sums(obj), 
                    estimate_richness(obj, measures="Shannon"),
                    btools::estimate_pd(obj),    #This adds Faith's phylo div and spec richness
                    "row"=rownames(sample_data(obj)),
                    evenness=(estimate_richness(obj, measures="Shannon")/log(estimate_richness(obj, measures = "Observed"))),
                    "cage_day"=paste0(sample_data(obj)$cage,sample_data(obj)$day),
                    "MT"= paste(sample_data(obj)$microb,sample_data(obj)$treatment, sep="_"))  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
exp_dat<- cbind(sample_data(psexp3000),dftmp)
head(exp_dat)

exp_dat$FinalReads <- exp_dat1$FinalReads #This is from me concatenating read counts outside of R. Pulling this into the workspace. 
exp_dat$MDT <- paste(sample_data(pstrans1)$div_day, sample_data(pstrans1)$treatment, sep="_")
exp_dat$microb <- factor(exp_dat$microb, levels=c("LO", "HI"))
exp_dat$treatment <- factor(exp_dat$treatment, levels=c("control", "worm"))
exp_dat$div_day <- factor(exp_dat$div_day, levels=c('LO_D-5', 'LO_D1', 'LO_D10', 'LO_D14', 'LO_D28', 
                                                    'LO_D31', 'LO_D41', 'HI_D-5', 'HI_D1', 'HI_D10', 
                                                    'HI_D14', 'HI_D28', 'HI_D31', 'HI_D41'))
exp_dat$MDT <- factor(exp_dat$MDT, 
                      levels=c('LO_D-5_worm', 'LO_D1_worm', 'LO_D10_worm', 'LO_D14_worm', 'LO_D28_worm','LO_D31_worm', 'LO_D41_worm',
                               'HI_D-5_worm', 'HI_D1_worm', 'HI_D10_worm', 'HI_D14_worm', 'HI_D28_worm', 'HI_D31_worm', 'HI_D41_worm', 
                               'LO_D-5_control', 'LO_D1_control', 'LO_D10_control', 'LO_D14_control', 'LO_D28_control', 'LO_D31_control', 'LO_D41_control', 
                               'HI_D-5_control', 'HI_D1_control', 'HI_D10_control', 'HI_D14_control', 'HI_D28_control', 'HI_D31_control', 'HI_D41_control'))
#Colors to use to define all control samples as grey, with variable being "MDT". 
allcols1 <- c("#e0c0f0",  '#d6a3f0' , '#c275eb' , '#A11FE5' , '#7E05BE' , '#5C008D', '#300049',
              "#faf1d7", "#f5e1a6" , "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00", 
              "#FFFFFF", "#FFFFFF", "#FFFFFF",  "#FFFFFF", "#FFFFFF", "#FFFFFF","#FFFFFF", "#FFFFFF", 
              "#FFFFFF", "#FFFFFF", "#FFFFFF",  "#FFFFFF", "#FFFFFF", "#FFFFFF")
#

#Subsetting df by treatment
treat_dat <- exp_dat[exp_dat$treatment=="worm",]
ctrl_dat <- exp_dat[exp_dat$treatment=="control",]


#Create a scatterplot of the richness and evenness data for all samples. 
g <- ggplot(exp_dat, aes(x=SR, y=evenness)) +
  geom_point(aes(fill=MDT, shape=microb), size=3) +
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values=allcols1) +
  #scale_colour_manual(values=allcols) +
  #scale_colour_manual(values=c('#62148c', '#b58e14')) +
  #scale_colour_manual(values=c('#FFFFFF', '#62148c', '#FFFFFF', '#b58e14')) +
  ggtitle("Richness vs. Evenness, all samples") +
  theme_bw() +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")


ggsave("AlphaDiv/psexp3000_RichnessEvennessScatter_allcolswhite.pdf", g, units="in", width=7, height=5)

#

# Grab the mins and maxes and sd's of each of the subgroups. 
obj <- ctrl_dat

df <- data.frame(minrichlo = min(obj$SR[obj$microb=="LO"]), 
                 maxrichlo = max(obj$SR[obj$microb=="LO"]), 
                 sdrichlo = sd(obj$SR[obj$microb=="LO"]),
                 minrichhi = min(obj$SR[obj$microb=="HI"]), 
                 maxrichhi = max(obj$SR[obj$microb=="HI"]), 
                 sdrichhi = sd(obj$SR[obj$microb=="HI"]))

df

#


#   ##### Generating distance matrices and plotting PCOAs #####

obj <- psexp3000
dftmp <- data.frame(estimate_richness(obj, measures = "Observed"),
                    sum = sample_sums(obj), 
                    estimate_richness(obj, measures="Shannon"),
                    btools::estimate_pd(obj),    #This adds Faith's phylo div and spec richness
                    "row"=rownames(sample_data(obj)),
                    evenness=(estimate_richness(obj, measures="Shannon")/log(estimate_richness(obj, measures = "Observed"))),
                    "cage_day"=paste0(sample_data(obj)$cage,sample_data(obj)$day),
                    "MT"= paste(sample_data(obj)$microb,sample_data(obj)$treatment, sep="_"))  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
exp_dat<- cbind(sample_data(psexp3000),dftmp)
head(exp_dat)

exp_dat$FinalReads <- exp_dat1$FinalReads #This is from me concatenating read counts outside of R. Pulling this into the workspace. 
exp_dat$MDT <- paste(sample_data(pstrans1)$div_day, sample_data(pstrans1)$treatment, sep="_")
exp_dat$microb <- factor(exp_dat$microb, levels=c("LO", "HI"))
exp_dat$treatment <- factor(exp_dat$treatment, levels=c("control", "worm"))
exp_dat$div_day <- factor(exp_dat$div_day, levels=c('LO_D-5', 'LO_D1', 'LO_D10', 'LO_D14', 'LO_D28', 
                                                    'LO_D31', 'LO_D41', 'HI_D-5', 'HI_D1', 'HI_D10', 
                                                    'HI_D14', 'HI_D28', 'HI_D31', 'HI_D41'))
exp_dat$MDT <- factor(exp_dat$MDT, 
                      levels=c('LO_D-5_worm', 'LO_D1_worm', 'LO_D10_worm', 'LO_D14_worm', 'LO_D28_worm','LO_D31_worm', 'LO_D41_worm',
                               'HI_D-5_worm', 'HI_D1_worm', 'HI_D10_worm', 'HI_D14_worm', 'HI_D28_worm', 'HI_D31_worm', 'HI_D41_worm', 
                               'LO_D-5_control', 'LO_D1_control', 'LO_D10_control', 'LO_D14_control', 'LO_D28_control', 'LO_D31_control', 'LO_D41_control', 
                               'HI_D-5_control', 'HI_D1_control', 'HI_D10_control', 'HI_D14_control', 'HI_D28_control', 'HI_D31_control', 'HI_D41_control'))
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
#Alternative D-5 color: #e0c0f0

#Make phylsoeq object compositional. 
pstrans1 <- transform_sample_counts(psexp3000, function(otu) otu/sum(otu))  

#Create ordinations with phyloseq. 
braypcoa <- ordinate(pstrans1, "PCoA", "bray")
unifpcoa <- ordinate(pstrans1, "PCoA", "unifrac")
wunifpcoa <- ordinate(pstrans1, "PCoA", "wunifrac")

#Generate GUniFrac distance matrix since you can't do this is phyloseq. 
obj <- psexp3000
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0, 0.5, 1))$unifracs #Generate the distances for all the weightings you want. 
d50 <- gunif[, , "d_0.5"] 

GU50pcoa <- data.frame(cmdscale(d50, eig=FALSE))
GU50pcoa <- GU50pcoa %>% 
  dplyr::rename("GU50_pc1" = "X1",
                "GU50_pc2" = "X2")
GU50pcoa <- GU50pcoa*(-1)
full <- cbind(exp_dat, GU50pcoa)

#

#    Plot the phyloseq ordinations # 

#Bray PCOA as the example. 
brayp <- plot_ordination(pstrans1, braypcoa) +
  geom_point(aes(colour=microb, fill=MDT, shape=treatment), size=3) +
  scale_shape_manual(values=c(21, 21)) +
  scale_fill_manual(values=allcols1) +
  scale_colour_manual(values=c('#7E05BE', '#b58e14')) +
  ylim(-0.45, 0.4)+
  xlim(-0.63, 0.4)+
  stat_ellipse(aes(color= microb), type="t", level=0.90)+
  ggtitle("Bray PCOA All Samples") +
  theme_bw() +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
brayp

#Also do this with the unweighted and weighted unifrac distances I generated in Phyloseq. 



#    GUniFrac 50% weighting PCOAs #
pc <- pcoa(d50) #generate the pcoa from your distance matrix so that you can grab the axis values to label the axes. 
GU50plot <- ggplot(aes(x=GU50_pc1, y=GU50_pc2, colour=label), data=full) +
  geom_point(aes(colour=microb, fill=MDT, shape=microb), size=3) +
  scale_shape_manual(values=c(21, 22)) +
  scale_fill_manual(values=allcols1) +
  #scale_colour_manual(values=c("#62148c","#62148c",'#b58e14', "#b58e14")) +
  scale_colour_manual(values=c('#62148c', '#b58e14')) +
  #scale_colour_manual(values=c('#FFFFFF', '#62148c', '#FFFFFF', '#b58e14')) +
  # ylim(-0.45, 0.4)+
  # xlim(-0.63, 0.4)+
  stat_ellipse(aes(color= microb), type="t", level=0.90)+
  theme_bw() +
  ggtitle("GUniFrac (50%) PCOA All samples") +
  xlab(paste("GUnif Axis 1 (", label_percent(accuracy=0.01)(pc$values[1,2]), ")", sep="")) + 
  ylab(paste("GUnif Axis 2 (", label_percent(accuracy=0.01)(pc$values[2,2]), ")", sep="")) +
  theme(plot.title=element_text(size=16, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        axis.text=element_text(size=15), 
        legend.position="none")
GU50plot

ggsave("BetaDiv/psexp3000treat_GUniFracPCOA_allcolswhite.pdf", GU50plot, units="in", width=7, height=5)




#Ok now regenerate the distance matrix for the worm-treated animals only. 
fulltreat <- full5[full5$treatment=="worm",]




############ Tracking the 11ish ASVs over time in both microbiome groups. #####


#   ##### Create the phyloseq object with the PhyLo11B taxa only. #####

#Outside of R, I looked at the ASVs that were found in psbubble and created a basic spreadsheet
# of the ASVs that were in at least half of samples. This totaled 15 ASVs that all identified
# to the 11 that we started with. there is a duplicate of Adlercreutzia and 3 ASVs of Muribaculum. 

#Read in those ASVs. 

orig <- read.csv("OTU_Tax_tables/Original_colonizers.csv")
#View(orig)
orig1 <- as.list(orig$ASV_ID)

#Ok now we want to subset our object of interest to only include those ASVs. 

tmp <- subset(otu_table(psexp), rownames(otu_table(psexp)) %in% orig1)
newobj <- merge_phyloseq(tmp, taxonomy=tax_table(psexp), metadata=sample_data(psexp), tree=phy_tree(psexp))
newobj #Seems like that worked! 15 ASVs and keeps all 324 samples. 

tmp <- subset(otu_table(psexp3000), rownames(otu_table(psexp3000)) %in% orig1)
newobj1 <- merge_phyloseq(tmp, taxonomy=tax_table(psexp3000), metadata=sample_data(psexp3000), tree=phy_tree(psexp3000))
newobj1

View(sample_data(newobj1))

#What do the read counts look like? Lower in the HI animals since they have so many more taxa. 
reads_unrarefied <- data.frame(readcount(newobj))

reads_unrarefied$AID <- rownames(reads_unrarefied)

reads_unrarefied <- reads_unrarefied %>%
  dplyr::rename(unrarefied = readcount.newobj.)

reads_rarefied <- data.frame(readcount(newobj1))

reads_rarefied$AID <- rownames(reads_rarefied)

reads_rarefied <- reads_rarefied %>% 
  dplyr::rename(rarefied = readcount.newobj1.)

new <- merge(reads_rarefied, reads_unrarefied, by="AID", all.y=TRUE)

View(new)
#Ok the rarefied and unrarefied give nearly the same order of representation of things, and 
# We should probably just stick with the rarefied since the unrarefied stuff seems to keep in
# some samples of very low quality. 


#Ok save the phyloseq objects. 
rare11 <- newobj1 #This is the rarefied data with only PhyLo11B taxa.
unrare11 <- newobj #This is the unrarefied data with only PhyLo11B taxa.


#   ##### Stacked bar of the PhyLo11B microbiomes in the bubble #####
psbubble1 <- prune_taxa(taxa_sums(psbubble)>13, psbubble)
psbubble1
psbubblefam <- tax_glom(psbubble1, taxrank="Family", NArm=TRUE)
phyloseq::taxa_names(psbubblefam) <- phyloseq::tax_table(psbubblefam)[, "Family"] #Change the names of the "ASVs" to the family name instead of representative QIIME2 ID. 
pseq <- microbiome::transform(psbubblefam, "compositional") #This is with all psprelo ASVs removed from all LO animals and psprehi removed from all HI animals. 

rownames(table(phyloseq::tax_table(pseq)[, "Family"]))

##Ok need to update colors to match the colors in the stacked area plot, if possible. 
p <- plot_composition(pseq, 
                      #group_by="day",
                      taxonomic.level="Family", 
                      otu.sort= "abundance",
                      transform = "compositional",
                      average_by = "day",
                      plot_type='barplot'
                      #sample.sort="sample.id",
) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=famcolsall) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=1),
        plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=14, face="bold"), 
        legend.text=element_text(size=13)) + 
  guides(fill=guide_legend(title="Family")) +
  labs(x = "Sampling Date", y = "Relative Abundance", title="LOW averaged rare Stacked bar, read counts", 
       subtitle="psbubble family stacked bar GG2 classification") + ggpubr::rotate_x_text(40) 
p

#
ggsave("Stacked_bars/psbubble_StackedBar_GG2classification_matchingColors.pdf", p, units="in", height=7, width = 4.5)
#



#   ##### Now let's look at a stacked bar to see what the relative representation looks like. #####


# Make sure we use functions from correct package
transform <- microbiome::transform


rare11@sam_data[["new"]] <- paste(rare11@sam_data[["microb"]], rare11@sam_data[["TD"]], sep="_")
rare11@sam_data[["MT"]] <- paste(rare11@sam_data[["microb"]], rare11@sam_data[["treatment"]], sep="_")
rare11@sam_data[["MT"]] <- factor(rare11@sam_data[["MT"]], levels=c("LO_control","LO_worm", "HI_control", "HI_worm"))

# Merge rare taxa to speed up examples
pseq <- transform(rare11, "compositional")
#Not going to aggregate for now because we want to see what each ASV is doing. 
pseq1 <- aggregate_rare(pseq, level = "Family", detection = 0.01, prevalence = 0.01)
pseq2 <- aggregate_rare(rare11, level = "Family", detection = 0.01, prevalence = 0.01)
#pseq1 <- top_taxa(pseq, n=12)

#Exploring the data
rownames(table(phyloseq::tax_table(psexp3000hi)[, "Family"]))
nrow(table(phyloseq::tax_table(pseq1)[, "Family"]))


pseq1@sam_data[["microb"]] <- factor(pseq1@sam_data[["microb"]], levels=c("LO","HI"))


#Make your colors work for you. 
nb.cols <- 11 #We know there are 15 ASVs here, so we can just keep them as such. 
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

p <- plot_composition(pseq2,
                      taxonomic.level="Family",
                      otu.sort= "abundance",
                      average_by = "new", 
                      group_by="MT",
                      transform = "compositional", 
                      #sample.sort="sample.id",
) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=mycolors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text.y=element_text(size=20)) +
  labs(x = "Sampling date", y = "Reads (out of 3000)", title="PhyLo11B ASVs averaged Stacked bar", 
       subtitle="1U4U") + ggpubr::rotate_x_text(40) 
p

print(p) +coord_flip() #If you want. 


#   ##### Pull the metadata. #####

dftmp <- data.frame(estimate_richness(rare11, measures = "Observed"),
                    estimate_richness(rare11, measures="Shannon"),
                    btools::estimate_pd(rare11),    #This adds Faith's phylo div and spec richness
                    "row"=rownames(sample_data(rare11)),
                    "evenness"=(estimate_richness(rare11, measures="Shannon")/log(estimate_richness(rare11, measures = "Observed"))),
                    "reads"=readcount(rare11),
                    #"MT"=paste0(sample_data(rare11)$treatment,sample_data(rare11)$microb) # Already in there from the last section. 
                    "cage_day"=paste0(sample_data(rare11)$cage,sample_data(rare11)$day) 
                    #"new" = paste(sample_data(rare11)$microb, sample_data(rare11)$TD, sep="_") # Already in there from the last section. 
)  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
dfrare11<- cbind(sample_data(rare11),dftmp)


#
#   ##### How many reads were retained just based on those ASVs throughout time? #####

#What do we want to do with our df?
View(dfrare11)
#Get the average number of reads and ASVs for each of the categories in the "new" column. 
# This separates both treatment and day. Could also just do the two columns adn get the average by
# aggregating. 

#THIS IS AN AMAZING FUNCTION THAT WOULD'VE SAVED ME SO MUCH TIME THROUGHOUT THE YEARS, WTF?!
avg1 <- aggregate(reads ~ microb + TD, data=dfrare11, FUN="mean")
avg <- aggregate(reads ~ new, data=dfrare11, FUN="mean")
View(avg1)
#Yield the same results, wow, this is amazing. 

mean(avg$reads[22:28])
sd(avg$reads[22:28])
#Hi_control: 201.54 (+/- 138.82) ... only 6.72% (+/- 4.63%) of total reads
#Hi_worm: 81.55 (+/- 16.30) ... 2.72% (+/- 0.54%) of total reads
#Low_worm: 2181.52 (+/- 205.36) ... 72.72% (+/- 6.85%) of total reads
#Low_control: 2055.04 (+/- 327.20) ... 68.50 (+/- 10.91%) of total reads


#

#      ##### Plot the number of reads represented just by the 15 ASVs. #####

dfrare11$microb <- factor(dfrare11$microb, level=c("LO", "HI"))

p1 <- ggplot(dfrare11, aes(x=day, y=reads)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) + 
  #geom_beeswarm(aes(shape=treatment, fill=microb), cex=1, size=8, alpha=0.9) +
  geom_point(aes(shape=treatment, fill=MT),
             position=position_jitterdodge(jitter.width=0.6), alpha=0.8, size=5) + 
  scale_colour_manual(values=c("#A11FE5", "#A11FE5","#DCAE1C", "#DCAE1C")) +
  scale_fill_manual(values=c("#A11FE5", "#A11FE5","#DCAE1C", "#DCAE1C")) +
  scale_shape_manual(values=c(24,21)) +
  labs(title= "Reads from PhyLo11B Synthetic Community", 
       subtitle = "",
       x="Sampling Day", y="Reads (out of 3000)") +
  theme_bw() +
  #facet_wrap(~microb) +
  theme(axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), legend.position="none")
p1




############ Remove all the HI taxa from the LOW dataset. #####

#Get the list of the pretreatment taxa from the HI animals. 
hitax <- as.list(rownames(otu_table(psexp3000hi)))

#Remove that from the LOW dataset on all days. 
tmpe <- subset(otu_table(psexp3000lo), !rownames(otu_table(psexp3000lo)) %in% hitax)
tmpe1 <- merge_phyloseq(tmpe, taxonomy=tax_table(psexp3000lo), metadata=sample_data(psexp3000lo))
tmpe1
#ok this still leaves us with 126 unique ASVs after psexp3000lo started with 
#This is enough to say that all these ASVs that are infiltrating are not coming
# from the HI animals. 


#      ##### Pull the metadata. #####

dftmp <- data.frame(estimate_richness(tmpe1, measures = "Observed"),
                    estimate_richness(tmpe1, measures="Shannon"),
                    # btools::estimate_pd(tmpe1),    #This adds Faith's phylo div and spec richness
                    "row"=rownames(sample_data(tmpe1)),
                    "evenness"=(estimate_richness(tmpe1, measures="Shannon")/log(estimate_richness(tmpe1, measures = "Observed"))),
                    "reads" = readcount(tmpe1), 
                    "cage_day"=paste0(sample_data(tmpe1)$cage,sample_data(tmpe1)$day), 
                    "treat_microb"=paste0(sample_data(tmpe1)$treatment,sample_data(tmpe1)$microb))  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
dfloNOhi<- cbind(sample_data(tmpe1),dftmp)

dfloNOhi$microb <- factor(dfloNOhi$microb, levels=c("LO", "HI"))
dfloNOhi$treatment <- factor(dfloNOhi$treatment, levels=c("worm", "control"))
dfloNOhi$TD <- factor(dfloNOhi$TD, levels=c('worm_D-5', 'worm_D1', 'worm_D10',  'worm_D14', 'worm_D28', 'worm_D31', 'worm_D41',
                                            'control_D-5', 'control_D1', 'control_D10',  'control_D14', 'control_D28', 'control_D31', 'control_D41'))


#      ##### How many of the ASVs are left? #####

p <- ggplot(dfloNOhi, aes(x=day_num, y=Observed)) +
  geom_boxplot(aes(color=TD), width=2.5, outlier.shape=NA) + 
  geom_beeswarm(aes(fill=TD), shape=21, dodge.width=2.5, alpha=0.8, size=5, na.rm = TRUE, cex=0.5) +  
  scale_fill_manual(values = locols1) +
  scale_color_manual(values = locolsx2) +
  #scale_shape_manual(values = c(21,21)) +
  theme_bw() +
  labs(title= "ASV richness in LOW but absent in HI", 
       subtitle = "Subtracting psexp3000hi ASVs from psexp3000lo",
       x="Sampling Day", y="ASV richness") +
  theme_bw() +
  theme(plot.title=element_text(size=20, face="bold"), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=18), 
        legend.position="none")
p
#Save it. 
#ggsave("AlphaDiv/psexp3000lo_WithoutHIASVs_RichnessBoxplots.pdf", units="in", height=7, width = 9)




############ Remove all pretreatment ASVs from each microbiome type and plot remainders #####

#Get the ASVs in the pretreatment samples for LO and HI.

lotax <- as.list(rownames(otu_table(psprelo)))
hitax <- as.list(rownames(otu_table(psprehi)))

#Create the subsetted phyloseq objects without D-5.

psa <- subset_samples(psexp3000lo, day != "D-5")
psa <- prune_taxa(taxa_sums(psa) >= 1, psa)

psb <- subset_samples(psexp3000hi, day != "D-5")
psb <- prune_taxa(taxa_sums(psb) >= 1, psb)


#Create the phyloseq objects that have all starting ASVs removed. 

length(lotax) #So this should remove 105 taxa from the updated phyloseq object. 
tmpa <- subset(otu_table(psa), !rownames(otu_table(psa)) %in% lotax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(psa), metadata=sample_data(psa))
#Actually removed 97 ASVs bc there were 8 ASVs that were in pretreatment but not in the other timepoints. 

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(psb), !rownames(otu_table(psb)) %in% hitax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psb), metadata=sample_data(psb))

psexp3000lo
psa
tmpa1


#Combine the LO and HI back together. 

psnopre1 <- merge_phyloseq(tmpa1, tmpb1, tree = phy_tree(psexp3000))
sample_data(psnopre1)$reads <- sample_sums(psnopre1)

#
#   ##### Pull the metadata and define factors. #####

dftmp <- data.frame(estimate_richness(psnopre1, measures = "Observed"),
                    estimate_richness(psnopre1, measures="Shannon"),
                    btools::estimate_pd(psnopre1),    #This adds Faith's phylo div and spec richness
                    "row"=rownames(sample_data(psnopre1)),
                    "evenness"=(estimate_richness(psnopre1, measures="Shannon")/log(estimate_richness(psnopre1, measures = "Observed"))),
                    "cage_day"=paste0(sample_data(psnopre1)$cage,sample_data(psnopre1)$day), 
                    "treat_microb"=paste0(sample_data(psnopre1)$treatment,sample_data(psnopre1)$microb))  
colnames(dftmp)[colnames(dftmp) == "Shannon.1"] ="evenness"
dfnopre1<- cbind(sample_data(psnopre1),dftmp)

dfnopre1$microb <- factor(dfnopre1$microb, levels=c("LO", "HI"))
dfnopre1$MT <- factor(paste(dfnopre1$microb, dfnopre1$treatment, sep="_"), levels=c("LO_control", "LO_worm", "HI_control", "HI_worm"))
dfnopre1$MTD <- paste(dfnopre1$microb, dfnopre1$day, dfnopre1$treatment, sep="_")
dfnopre1$MTD <- factor(dfnopre1$MTD, 
                       levels=c('LO_D1_worm', 'LO_D10_worm', 'LO_D14_worm', 'LO_D28_worm','LO_D31_worm', 'LO_D41_worm',
                                'HI_D1_worm', 'HI_D10_worm', 'HI_D14_worm', 'HI_D28_worm', 'HI_D31_worm', 'HI_D41_worm', 
                                'LO_D1_control', 'LO_D10_control', 'LO_D14_control', 'LO_D28_control', 'LO_D31_control', 'LO_D41_control', 
                                'HI_D1_control', 'HI_D10_control', 'HI_D14_control', 'HI_D28_control', 'HI_D31_control', 'HI_D41_control'))

dfnopre1treat <- dfnopre1[dfnopre1$treatment=="worm",]
dfnopre1ctrl <- dfnopre1[dfnopre1$treatment=="control",]

#
#      ##### Write out the OTU/Taxa table. #####
psnopre1_hictrl <- subset_samples(psnopre1, treatment=="control" & microb=="HI")
psnopre1_hictrl <- prune_taxa(taxa_sums(psnopre1_hictrl) > 1, psnopre1_hictrl)

#To pull out the taxa and otu tables and plot outside of phyloseq:
OTU <- data.frame(otu_table(psnopre1_hictrl))
OTU$sum <- rowSums(OTU)
otu_tax <- cbind(OTU, tax_table(psnopre1_hictrl), 
                 "ASV"=rownames(OTU))

#
#   ##### Plot reads, alpha diversity of the outcomes #####

#Plot the number of ASVs remaining after removing all pretreatmeent ASVs from LO and HI
p1 <- ggplot(dfnopre1, aes(x=day_num, y=Observed)) +
  geom_boxplot(aes(color=MTD), width=3, outlier.shape=NA) +
  geom_beeswarm(aes(fill=MTD, shape=microb), dodge.width=3, alpha=0.8, size=5, na.rm = TRUE, cex=1.1) +  
  scale_colour_manual(values=allcolsdiffgrey) +
  scale_fill_manual(values=allcolsdiff) +
  scale_shape_manual(values=c(21,22)) +
  labs(title= "ASVs added after pretreatment", 
       subtitle = "Subtracting same ASVs from both worm-treated and controls",
       x="Sampling Day", y="# of ASVs") +
  theme_bw() +
  facet_wrap(~microb) +
  theme(axis.title=element_text(size=20, face="bold"), 
        axis.text=element_text(size=20), legend.position="none")
p1


ggsave("AlphaDiv/psexp3000_nopretreatASVsBoxplots_allcolswhite.pdf", p1, units="in", width=14, height=7)


#       #Summarize the data #####
summary(dfnopre1ctrl$reads[dfnopre1ctrl$microb=="LO" & dfnopre1ctrl$day=="D10"])
sd(dfnopre1ctrl$reads[dfnopre1ctrl$microb=="LO" & dfnopre1ctrl$day=="D10"])
#
summary(dfnopre1treat$Observed[dfnopre1treat$microb=="LO" & dfnopre1treat$day=="D41"])
sd(dfnopre1treat$Observed[dfnopre1treat$microb=="LO" & dfnopre1treat$day=="D41"])

#


#   ##### Model these outcomes. #####

# Ok, now model the differences between timepoints and microbiome types. # 

mem.invaders <- lme(Observed ~ day_num * microb + treatment * day_num + microb * treatment, random=~1|mouse, data=dfnopre1)
summary(mem.invaders)


mem.worminvaders <- lme(Observed ~ day_num * microb, random=~1|mouse, data=dfnopre1treat)
summary(mem.worminvaders)

#   ##### Create a stacked bar plot of the invaders. #####

psnopre1@sam_data[["MT"]] <- paste(psnopre1@sam_data[["microb"]], psnopre1@sam_data[["treatment"]], sep="_")
psnopre1@sam_data[["MT"]] <- factor(psnopre1@sam_data[["MT"]], levels=c("LO_control","LO_worm", "HI_control", "HI_worm"))
psnopre1@sam_data[["new"]] <- paste(psnopre1@sam_data[["microb"]], psnopre1@sam_data[["TD"]], sep="_")

psnopre@sam_data[["MT"]] <- paste(psnopre@sam_data[["microb"]], psnopre@sam_data[["treatment"]], sep="_")
psnopre@sam_data[["MT"]] <- factor(psnopre@sam_data[["MT"]], levels=c("LO_control","LO_worm", "HI_control", "HI_worm"))
psnopre@sam_data[["new"]] <- paste(psnopre@sam_data[["microb"]], psnopre@sam_data[["TD"]], sep="_")


#Create phyloseq objects to look at worm and sham separately. 
psnopretreat <- subset_samples(psnopre1, treatment=="worm")
psnopretreat <- prune_taxa(taxa_sums(psnopretreat) >= 1, psnopretreat)

psnoprectrl <- subset_samples(psnopre1, treatment=="control")
psnoprectrl <- prune_taxa(taxa_sums(psnoprectrl) >= 1, psnoprectrl)

# Make sure we use functions from correct package
transform <- microbiome::transform


# Merge rare taxa to speed up examples
pseqlo <- transform(tmpa1, "compositional") #LO samples.
pseqhi <- transform(tmpb1, "compositional") #HI samples. 
#Not going to aggregate for now because we want to see what each ASV is doing. 
pseq <- microbiome::transform(psnopre1, "compositional") #This is with all psprelo ASVs removed from all LO animals and psprehi removed from all HI animals. 
pseq <- microbiome::transform(psnopre, "compositional") #This is with all pretreatment ASVs removed from each MT group (e.g., psprelotreat removed from LO treat)
pseq1 <- aggregate_rare(psnopretreat, level = "Family", detection = 0.01, prevalence = 0.05)
pseq2 <- aggregate_rare(pseqhi, level = "Family", detection = 0.01, prevalence = 0.05)


#Exploring the data
rownames(table(phyloseq::tax_table(pseq1)[, "Family"]))
sum(colSums(otu_table(pseq1))) #How many samples actually have some of these ASVs represented?
ncol(otu_table(pseq1)) #The difference between this and the last line will tell us how many
# samples have none of the common ASVs infiltrating.

#Organize the data for faceting, if need be. 
pseq1@sam_data[["microb"]] <- factor(pseq1@sam_data[["microb"]], levels=c("LO","HI"))
pseq2@sam_data[["treatment"]] <- factor(pseq2@sam_data[["treatment"]], levels=c("worm","control"))
pseq3@sam_data[["microb"]] <- factor(pseq3@sam_data[["microb"]], levels=c("LO","HI"))


p <- plot_composition(pseq1, 
                      group_by="MT",
                      #taxonomic.level="Family", 
                      #otu.sort= "abundance",
                      transform = "compositional",
                      average_by = "new",
                      plot_type='heatmap'
                      #sample.sort="sample.id",
) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=famcolsall) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=1),
        plot.title=element_text(size=16, face="bold"), 
        axis.title=element_text(size=15, face="bold"), 
        axis.text=element_text(size=15), 
        legend.title=element_text(size=14, face="bold"), 
        legend.text=element_text(size=13)) +   
  labs(x = "Microbiome Type / Timepoint", y = "Reads (out of 3000)", title="Averaged rare Stacked bar, read counts, Treatment Only", 
       fill= "Family",
       subtitle="After removing same pretreatment ASVs from Worm and Sham (detection = 0.01, prevalence = 0.05)") + ggpubr::rotate_x_text(40) 
p


#Save this plot. 
ggsave("Stacked_bars/psexp3000treat_NoPretreatASVsStackedBar_InvadingASVs_ReadCounts.pdf", p, units="in", height=7, width=11)

#Not all bars are going to be complete since not all animals have these infiltrators at all days. 
print(p) +coord_flip() #If you want. 


#
#   ##### Now, how many ASVs of invaders are there in each group? #####

avg <- aggregate(Observed ~ new, data=dfnopre1, FUN="mean")
View(avg)

mean(avg$Observed[1:6])
sd(avg$Observed[1:6])


#

############ Create phylogeny of all Lachnospiraceae ASVs in LO and HI #####
#   ##### Create an average, compositional phyloseq object #####
pscomp <- microbiome::transform(psexp3000treat, "compositional")
pslach <- subset_taxa(pscomp, Family == "Lachnospiraceae")

#Grab taxa table and tree. 
TAX <- phyloseq::tax_table(pslach)
TREE <- phy_tree(pslach)

#Now prep the otu_table. 
#Grab all your taxa from the hi animals, average them. 
pslach.hi <- subset_samples(pslach, microb=="HI")
#pslach.hi1 <- subset_samples(pslach, mouse %in% c("14_RLP", "17_LLP"))
lachavg.hi <- data.frame(rowMeans(otu_table(pslach.hi)))
lachavg.hi <- lachavg.hi %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    HI=rowMeans.otu_table.pslach.hi..)

#Grab all your taxa from the LO animals, average them. 
pslach.lo <- subset_samples(pslach, microb=="LO")
lachavg.lo <- data.frame(rowMeans(otu_table(pslach.lo)))
lachavg.lo <- lachavg.lo %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    LO=rowMeans.otu_table.pslach.lo..)
lachavg <- cbind(lachavg.hi, lachavg.lo)

#Make it officially an OTU table. 
OTU <- otu_table(lachavg, taxa_are_rows = TRUE)

#Make your teeny tiny metadata with whatever variables you want to use to label tree. 
microb=c("LO", "HI")
meta <- data.frame(microb)
rownames(meta) <- c("LO", "HI")

#Merge it all back into a phyloseq object. 
pslach.avg <- phyloseq(OTU, TAX, TREE)
pslach.avg1 <- merge_phyloseq(pslach.avg, sample_data(meta))

#   ##### Create a phyloseq object that is subsetted to one family and based on metadata #####
pslach.hi <- subset_samples(pslach, microb=="HI")
pslach.lo <- subset_samples(pslach, microb=="LO")

#Now remerge after doing whatever you needed to do.  
pstmp <- merge_phyloseq(pslach.lo, pslach.hi)
pstmp <- prune_taxa(taxa_sums(pstmp)>1, pstmp)


#       ##### Make the tree #####

lach.plot <- plot_tree(pslach.avg1, color="microb", label.tips = "Genus", shape="microb",
                       size = "abundance", plot.margin = 0.3, text.size=3, ladderize = TRUE) +
  #scale_color_brewer(palette="Blues") +
  scale_shape_manual(values=c(22, 21)) +
  scale_fill_manual(values = c("#DCAE1C",'#A11FE5')) +
  scale_color_manual(values=c("black", "black")) +
  ggtitle("psexp3000treat, Lachnospiraceae average relative abundances")

lach.plot




############ Create a phylogeny of the genera that were present in the pretreatment and posttreatment samples. #####
# This is specifically the code for only treatment or only control samples #
#Phylogenetic Network Plot for only treatment animals. 
#Starting with the same workspace setup as the main Phylogenetic Network plot script. 

#https://github.com/lch14forever/microbiomeViz/tree/master

#For the phylogeny of all taxa in the group, and then those that are marked based on day and group. 

#
psexp3000

#Define important commands that have been overwritten by other packages. 
rename <- dplyr::rename
tax_table <- phyloseq::tax_table


#Subset to just the group you want. 
pstmp <- subset_samples(psexp3000, treatment=="worm") 
psexp3000treat <- prune_taxa(taxa_sums(pstmp)>1, pstmp)

pstmp <- subset_samples(psexp3000, treatment=="control") 
psexp3000ctrl <- prune_taxa(taxa_sums(pstmp)>1, pstmp)


#Or can condense to genus or family level. 
#pssubfam <- tax_glom(pssub, taxrank=rank_names(pssub)[5]) #This was just a couple cages to test the code. 
psfam <- tax_glom(psexp3000treat, taxrank=rank_names(psexp3000treat)[5]) 
psgen <- tax_glom(psexp3000treat, taxrank=rank_names(psexp3000ctrl)[6]) #Use the whole phyloseq object to create the backbone. 

phyloseq::taxa_names(psfam) <- phyloseq::tax_table(psfam)[, "Family"] #Change the names of the "ASVs" to the family name instead of representative QIIME2 ID. 
phyloseq::taxa_names(psgen) <- phyloseq::tax_table(psgen)[, "Genus"]


#Working from microbiomeViz start and then with Zhou GitHub code. 
pstrans1 <- transform_sample_counts(psgen, function(otu) otu/sum(otu))
pstrans1 = filter_taxa(pstrans1, function(x) max(x)>=0.001,TRUE)
pstrans1 = fix_duplicate_tax(pstrans1)

tr = parsePhyloseq(pstrans1, use_abundance=FALSE, node.size.scale=0, node.size.offset=0) #Using microbiomeViz package for these. 
raw_p = tree.backbone(tr, size=0.5, fill="black", shape=16, color='black', layout="circular") + 
  ggtitle("Tree backbone —- Genus level") #Ok this makes a basic tree. Need to now add annotations to it. 
raw_p$data$nodeSize = 0.05
raw_p

#Define the variable that we're going to use to change node size on the plot.  
raw_p$data$nodeClass2 = as.character(raw_p$data$nodeClass)
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "f"] = "Family"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "c"] = "Class"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "o"] = "Order"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "p"] = "Phylum"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "g"] = "Genus"
raw_p$data$nodeClass2[raw_p$data$nodeClass2 == "k"] = "Kingdom"

#Changing the size of the circles around the nodes. 
raw_p$data$nodeSize = 0.2
raw_p$data$nodeSize2 = raw_p$data$nodeSize
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Root"] = 2
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Kingdom"] = 1.5
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Phylum"] = 1.25
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Class"] = 1
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Order"] = 0.75
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Family"] = 0.50
raw_p$data$nodeSize2[raw_p$data$nodeClass2 == "Genus"] = 0.25

#Ok now update the tree. 

#If you want to have color on the internal nodes:
#raw_p =
#  raw_p +
#  geom_point2(aes(color = nodeClass2,
#                  size = nodeSize2),
#              show.legend = FALSE) +
#  ggsci::scale_color_jama() +
#  ggnewscale::new_scale(new_aes = "fill")
#raw_p

#Or just make it with black internal nodes. 
raw_p =
  raw_p +
  geom_point2(aes(size = nodeSize2), color="black",
              show.legend = TRUE) 
#  ggsci::scale_color_jama() 
#  ggnewscale::new_scale(new_aes = "fill")

raw_p #This is your tree. You can look at all the layers with raw_p$ and then any of the info that pops up.  

#Pulling the names of the taxa from the tree. 
raw_p$data$label2 =
  raw_p$data$label %>%
  stringr::str_split("__") %>%
  purrr::map(function(x) {
    x[2]
  }) %>%
  unlist()

#Remove all names that aren't the level of focus. 
raw_p$data$label2[as.character(raw_p$data$nodeClass) != "g"] = NA  #This would be "f" for family. 

raw_p$data$label2
#raw_p$data
#



#   ##### If you're assembling from multiple phyloseq objects, this is how to pull each apart and then join them #####
# with only the unique tips you're interested in. ##
#Ok, now we need to curate and add the other data that we're going to put into the heatmap. 
#Essentially need to make a list of the families present in each of the important groups, so probably
# pretreatment LO and HI, and Cleared LO and HI. 
#First subset the phyloseq objects
#These all say "fam" because I started with family level. However, doing other levels too. Just leaving the labels to make things easy. 

pstmp <- subset_samples(pstrans1, microb=="LO" & day=="D-5" & treatment =="worm") #Subset and trim.
pssubfam1 <- prune_taxa(taxa_sums(pstmp)>0, pstmp)
loearly_expression_dat <- pssubfam1@otu_table@.Data %>% as.data.frame() #Expression data means otu_table, which is OTUs and read counts
loearly_variable_info <- as.data.frame(pssubfam1@tax_table@.Data) #Variable info means tax_table, which is all taxa names. 
loearly_sample_info = get_variable(physeq = pssubfam1)#Sample info is sample_data, which is all metadata. 

pstmp <- subset_samples(pstrans1, microb=="LO" & day=="D41" & treatment =="worm") #Subset and trim.
pssubfam2 <- prune_taxa(taxa_sums(pstmp)>0, pstmp)
lolate_expression_dat <- pssubfam2@otu_table@.Data %>% as.data.frame()
lolate_variable_info <- as.data.frame(pssubfam2@tax_table@.Data) 
lolate_sample_info = get_variable(physeq = pssubfam2)

pstmp <- subset_samples(pstrans1, microb=="HI" & day=="D-5" & treatment =="worm") #Subset and trim.
pssubfam3 <- prune_taxa(taxa_sums(pstmp)>0, pstmp)
hiearly_expression_dat <- pssubfam3@otu_table@.Data %>% as.data.frame()
hiearly_variable_info <- as.data.frame(pssubfam3@tax_table@.Data) 
hiearly_sample_info = get_variable(physeq = pssubfam3)

pstmp <- subset_samples(pstrans1, microb=="HI" & day=="D41" & treatment =="worm") #Subset and trim.
pssubfam4 <- prune_taxa(taxa_sums(pstmp)>0, pstmp)
hilate_expression_dat <- pssubfam4@otu_table@.Data %>% as.data.frame()
hilate_variable_info <- as.data.frame(pssubfam4@tax_table@.Data) 
hilate_sample_info = get_variable(physeq = pssubfam4)

#   ##### Put all the variable info (i.e., taxa tables) together. #####
#This just pulls all the taxa tables back together after pulling them apart. 
#Since I started with it all together, this is the same as pulling the taxa table from the original pssubfam. 
variable_info =
  dplyr::full_join(loearly_variable_info,
                   lolate_variable_info,
                   by = colnames(loearly_variable_info)) %>%
  dplyr::full_join(hiearly_variable_info,
                   by = colnames(loearly_variable_info)) %>%
  dplyr::full_join(hilate_variable_info,
                   by = colnames(loearly_variable_info)) %>%
  dplyr::distinct(Genus, .keep_all = TRUE)
#
#   ##### Didn't end up using this code, but he also created dfs to create a combined physeq object. #####
# Rename the rows to have ASV numbers. 
rownames(variable_info) = paste("asv", 1:nrow(variable_info), sep="_")
rownames(variable_info)[1:5]


# This just creates a shell matrix the the same number of rows as the number of families. 
#I didn't use this chunk that's commented out...
#expression_data =
#  matrix(1:(nrow(variable_info) * 2), ncol = 2)

#colnames(expression_data) = paste("sample", 1:ncol(expression_data), sep = "_")
#rownames(expression_data) = rownames(variable_info)

#sample_info =
#  data.frame(sample_id = colnames(expression_data))
#rownames(sample_info) = sample_info$sample_id

#expression_data = otu_table(expression_data, taxa_are_rows = TRUE)
#variable_info = tax_table(as.matrix(variable_info))
#sample_info = sample_data(sample_info)

#physeqGenus = phyloseq(expression_data, variable_info, sample_info)

#physeqGenus2 =
#  subset_taxa(physeq = physeqGenus, Genus %in% temp_data$genus) #temp_data is just the 
#df that had all the unique genera  
#and dissimilarity scores in the original Zhou et al. paper. 



#
#   ##### Creating temp df to include the annotations for each of the genera. #####
#So I should create a df with the family names and some score of whether they have the family or not. 


#So my rows for my version of "stool_personalized_score" could be something like: family, score (binary 0 or 1), dataset 
# (the last being the group category, like "lolate"). Have one for each physeq subobject. Then bind them into temp_data
loearly_tmp <- data.frame("Genus" = loearly_variable_info$Genus, 
                          "score" = 1, 
                          "class" = "loearly") #Try something like this for the other 3 subobjects and then move to temp_data creation and forward...
lolate_tmp <- data.frame("Genus" = lolate_variable_info$Genus, 
                         "score" = 1, 
                         "class" = "lolate")
hiearly_tmp <- data.frame("Genus" = hiearly_variable_info$Genus, 
                          "score" = 1, 
                          "class" = "hiearly")
hilate_tmp <- data.frame("Genus" = hilate_variable_info$Genus, 
                         "score" = 1, 
                         "class" = "hilate")

temp_data =  #Bind all the tmp df's together. 
  rbind(
    loearly_tmp,
    lolate_tmp,
    hiearly_tmp,
    hilate_tmp
  ) %>%
  dplyr::mutate(class = stringr::str_to_sentence(class)) #Last command just changes the first letter of the first word to uppercase. 


#   ##### Create a Venn Diagram showing overlapping ASVs (or families in this case) ########
mycols <- c("#e0c0f0",  '#5C008D', "#faf1d7",  "#8F6D00" )

Vennplot =
  ggVennDiagram(
    x =
      list(
        loearly = unique(loearly_tmp$Genus),
        lolate = unique(lolate_tmp$Genus),
        hiearly = unique(hiearly_tmp$Genus),
        hilate = unique(hilate_tmp$Genus)
      ),
    label_color = "white",
    label_geom = "text"
  ) +
  scale_colour_manual(values = c("#e0c0f0",  '#5C008D', "#faf1d7",  "#8F6D00" )) #Colors aren't working here, but probably not coloring right. Seems like it should be a gradient. 
Vennplot


#
#   ##### Now prep the data to populate the beautiful phylogeny. #####
#Ok unmelt the df so the it's wide instead of long. Doing this with tidyr method instead of reshape2. 
temp_data_scores = temp_data %>%
  dplyr::select(Genus, class, score) %>%
  tidyr::pivot_wider(names_from = class, values_from = "score")

temp_data_scores1 <- temp_data_scores %>%   #Duplicate this so that you can have another one with 
  rename(Loearly_score = Loearly, 
         Lolate_score = Lolate, 
         Hiearly_score = Hiearly, 
         Hilate_score = Hilate) 
temp_data_scores1[is.na(temp_data_scores1)] = 0

temp_data_scores1
#
#   ##### Ok now to add the metadata into a df that we can combine with the tree. ######
#At line 390 on the Zhou et al. 2024 GitHub page. 
#Make sure to update the level of taxonomy if not using Genus level. 
new_info =
  data.frame(Genus = raw_p$data$label2) %>%
  dplyr::left_join(as.data.frame(variable_info)[, c("Genus", "Kingdom", "Phylum", "Class", "Order", "Family")], #This is the same as calling the tax_table from the whole phyloseq object since I'm starting with one whole object. 
                   by = "Genus") %>%
  dplyr::left_join(temp_data_scores, by = c("Genus" = "Genus")) %>% #I don't think I need to do this...This is to add multiple annotating factors into the same df. The only other thing I could see doing is to pull the LFC or p-adj from DESeq2 and annotating the families that include significant ASVs enriched in these groups. 
  dplyr::mutate(Loearly = case_when(!is.na(Loearly) ~ "LO Early", TRUE ~ "no")) %>%
  dplyr::mutate(Lolate = case_when(!is.na(Lolate) ~ "LO Late", TRUE ~ "no")) %>%
  dplyr::mutate(Hiearly = case_when(!is.na(Hiearly) ~ "HI Early", TRUE ~ "no")) %>%
  dplyr::mutate(Hilate = case_when(!is.na(Hilate) ~ "HI Late", TRUE ~ "no")) %>%
  dplyr::left_join(temp_data_scores1, by = c("Genus" = "Genus"))

new_info
#
#         

#   ##### Pull it together with the tree data and define your colors.  #####
raw_p$data = cbind(raw_p$data, new_info)
#

#add tip label (genus label) and add point to show this body is here or not
mycols <- c("#e0c0f0",  '#5C008D', "#faf1d7",  "#8F6D00" )

mycols1 = c( "`Lo Early`" = "#e0c0f0",
             "`Lo Late`" = '#5C008D',
             "`Hi Early`" = "#faf1d7",
             "`Hi Late`" = "#8F6D00")

#   ##### Make a tree with just presence/absence noted on the tips with grey/white circles. #####

#This isn't super necessary, but it does give a presence/absence of the genera in each group. 
# Following up with heatmaps that are more useful in a couple more steps. 
#Colors aren't working, but the "no" does come up as white. Interesting if just looking at presence/absence. 
raw_p1 =
  raw_p +
  geom_tippoint(
    aes(fill = Loearly),
    shape = 21,
    x = 6.2,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(mycols1, "no" = "white")) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    aes(fill = Lolate),
    shape = 21,
    x = 6.4,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(mycols1[2], "no" = "white")) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    aes(fill = Hiearly),
    shape = 21,
    x = 6.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(mycols1, "no" = "white")) +
  ggnewscale::new_scale(new_aes = "fill") +
  geom_tippoint(
    aes(fill = Hilate),
    shape = 21,
    x = 6.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(mycols1, "no" = "white")) +
  ggnewscale::new_scale(new_aes = "fill")

raw_p1


#Copy to mess around with colors
raw_p1 =
  raw_p +
  geom_tippoint(
    aes(fill = Loearly),
    shape = 21,
    x = 6.2,
    size = 3,
    show.legend = FALSE
  ) +
  geom_tippoint(
    aes(fill = Lolate),
    shape = 21,
    x = 6.4,
    size = 3,
    show.legend = FALSE
  ) +
  geom_tippoint(
    aes(fill = Hiearly),
    shape = 21,
    x = 6.6,
    size = 3,
    show.legend = FALSE
  ) +
  geom_tippoint(
    aes(fill = Hilate),
    shape = 21,
    x = 6.8,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(mycols1, "no" = "white")) +
  ggnewscale::new_scale(new_aes = "fill") 

raw_p1

#   ##### add highlight for the phylum #####
phylum_color <- colorRampPalette(brewer.pal(8, "Set2"))(10) #Grabs 10 colors from across the whole Set2 palette
nrow(tax_table(pstrans1)[2]) 

hight_data =
  raw_p1$data %>%
  dplyr::filter(stringr::str_detect(label, "p__")) %>%
  dplyr::select(node, label) %>%
  dplyr::mutate(label = stringr::str_replace_all(label, "p__", "")) %>%
  dplyr::rename(id = node, type = label)

raw_p2 = raw_p +
  geom_hilight(data = hight_data, mapping = aes(node = id, fill = type), alpha = .4) +
  scale_fill_manual(values=phylum_color)

raw_p2 #This has the phylum annotations as well as whether the family is present in the sample. 


#Ok now need to add the heatmaps for 0 and 1 and the taxa names. 


##add multiple tip information
Loearly_score_info =
  raw_p2$data[, c("Loearly_score", "isTip")]
rownames(Loearly_score_info) = raw_p2$data$label
Loearly_score_info = Loearly_score_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

Lolate_score_info =
  raw_p2$data[, c("Lolate_score", "isTip")]
rownames(Lolate_score_info) = raw_p2$data$label
Lolate_score_info = Lolate_score_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

Hiearly_score_info =
  raw_p2$data[, c("Hiearly_score", "isTip")]
rownames(Hiearly_score_info) = raw_p2$data$label
Hiearly_score_info = Hiearly_score_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)

Hilate_score_info =
  raw_p2$data[, c("Hilate_score", "isTip")]
rownames(Hilate_score_info) = raw_p2$data$label
Hilate_score_info = Hilate_score_info %>%
  dplyr::filter(isTip) %>%
  dplyr::select(-isTip)


range(
  c(
    raw_p2$data$Loearly_score,
    raw_p2$data$Lolate_score,
    raw_p2$data$Hiearly_score,
    raw_p2$data$Hilate_score
  ),
  na.rm = TRUE
)

p1 = raw_p2 +
  ggnewscale::new_scale_fill()


#
p2_1 =
  gheatmap(
    p = p1,
    data = Loearly_score_info,
    offset = 0.05,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = "black",
    legend_title = "Index1"
  ) +
  scale_fill_gradient(
    low= "white",
    high= mycols[1],
    na.value = "white",
    limits = c(0, 1)
  ) +
  ggnewscale::new_scale(new_aes = "fill") 
#  geom_tippoint(
#    mapping = aes(shape = Stool_fc1_star),
#    x = 6.35,
#    size = 2,
#    show.legend = FALSE
#  ) +
#  scale_shape_manual(values = c("*" = "*")) +
#  ggnewscale::new_scale(new_aes = "shape")

p2_1
#



#
p2_2 =
  gheatmap(
    p = p2_1,
    data = Lolate_score_info,
    offset = 0.55,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = "black",
    legend_title = "Index1"
  ) +
  scale_fill_gradient(
    low= "white",
    high= mycols[2],
    na.value = "white",
    limits = c(0, 1)
  ) +
  ggnewscale::new_scale(new_aes = "fill") 
p2_2
#



#
p2_3 =
  gheatmap(
    p = p2_2,
    data = Hiearly_score_info,
    offset = 1.05,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = "black",
    legend_title = "Index1"
  ) +
  scale_fill_gradient(
    low= "white",
    high= mycols[3],
    na.value = "white",
    limits = c(0, 1)
  ) +
  ggnewscale::new_scale(new_aes = "fill") 
p2_3





p2_4 =
  gheatmap(
    p = p2_3,
    data = Hilate_score_info,
    offset = 1.55,
    width = .08,
    colnames_angle = 95,
    colnames_offset_y = .25,
    colnames = FALSE,
    color = "black",
    legend_title = "Index1"
  ) +
  scale_fill_gradient(
    low= "white",
    high= mycols[4],
    na.value = "white",
    limits = c(0, 1)
  ) +
  ggnewscale::new_scale(new_aes = "fill") 

p2_4





##add tip lab
##only add some tip points
# idx1 = which(!is.na(raw_p$data$Stool_fc1_p_adjust) & raw_p$data$Stool_fc1_p_adjust < 0.001)
# idx1 = sample(idx1, 10)
# raw_p$data$label2[-idx1] = NA
p3 =
  p2_4 +
  ggnewscale::new_scale(new_aes = "color") +
  geom_tiplab(
    aes(label = label2,
        color = Phylum),
    offset = 2.5,
    size =5,
    show.legend = FALSE
  ) +
  scale_color_manual(values=phylum_color) +
  ggtitle("Genus level Tree with labels <-  Control animals only") +
  theme(legend.position="none")

p3


#ggsave(p3, filename = "psexp3000treat_draft_cladogram_genus_nolegend.pdf", width = 12, height = 10, scale=1, limitsize=FALSE)
ggsave(p3, filename = "psexp3000ctrl_draft_cladogram_genus_nolegend.pdf", width = 12, height = 10, scale=1, limitsize=FALSE)



#


#


#         ##### IMPORTANT PLOTS: Bray-Curtis #####

mem.hitreat <- lmer(bray~ day_num_1 + (1|mouse_1), data=dist.dm.time.hi)
mem.lotreat <- lmer(bray~ day_num_1 + (1|mouse_1), data=dist.dm.time.lo)
mem.losham <- lmer(bray~ day_num_1 + (1|mouse_1), data=losham)
mem.hisham <- lmer(bray~ day_num_1 + (1|mouse_1), data=hisham)

fixef(mem.lotreat)
fixef(mem.hitreat)
fixef(mem.losham)
fixef(mem.hisham)


p12 = ggplot(dist.dm.time, aes(x=day_num_1, y=bray)) + 
  #geom_boxplot(aes(color=MDT), width=4, outlier.shape=NA) +
  geom_beeswarm(aes(fill=MDT, shape=microb_1), dodge.width=3, alpha=0.8, size=5, na.rm = TRUE, cex=1.1) +  
  #geom_point(aes(shape=treatment_1, fill=MDT), size=2, alpha=0.8, position=position_dodge(width=4)) +
  #theme_biome_utils() + 
  ylab("Individual dissimilarity score") + 
  xlab("Days post-infection") + 
  scale_shape_manual(values=c(21,22)) + 
  scale_fill_manual(values=allcols1) +
  #scale_fill_manual(values=c('#6C696E', '#A11FE5','#6C696E', "#DCAE1C")) + 
  #scale_color_manual(values=allcols1) +
  # geom_smooth(method="lm", aes(x=day_num_1, y=unifrac, color=microb_1), level=0.9, alpha=0.2) + 
  #scale_color_manual(values=c('#A11FE5',"#DCAE1C")) +
  #scale_color_identity() +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size=20), 
        axis.title=element_text(size=20, face="bold"), 
        axis.text.y=element_text(size=20)) +
  #stat_poly_eq(use_label(c("eq", "R2"))) +  
  facet_wrap(~microb_1) +
  ggtitle("Bray-Curtis Change since pretreat")
#ok now make sure you have your lmer results. Use "Estimate" column.
p12allcols <- p12 + 
  geom_abline(intercept=data.frame(fixef(mem.lotreat))[1,1], slope=data.frame(fixef(mem.lotreat))[2,1], colour="#A11FE5") + #Outcome of lotreat only model
  geom_abline(intercept=data.frame(fixef(mem.hitreat))[1,1], slope=data.frame(fixef(mem.hitreat))[2,1], colour="#DCAE1C") + #Outcome of hitreat
  geom_abline(intercept=data.frame(fixef(mem.losham))[1,1], slope=data.frame(fixef(mem.losham))[2,1], colour="#383838") + #Outcome of losham
  geom_abline(intercept=data.frame(fixef(mem.hisham))[1,1], slope=data.frame(fixef(mem.hisham))[2,1], colour="#292929")  #Outcome of hisham
p12allcols




