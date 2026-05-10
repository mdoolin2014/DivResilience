# Script to run richness and beta diversity change.
# Separated from the longer script for the whole manuscript for ease of tracking work. 
#Duplicated from original on 18Nov2025 because we think we finally have decided on the plots and 
# analyses we want to publish. 
#By Maggie Doolin. Has been cobbled from multiple sources. 

##### Set up workspace, load libraries, define colors #####

setwd("~/Documents/1U4U_16S/GG2_reclassification")
#load packages. 
library(ggplot2)
library(ggbeeswarm)
library(ape)
library(dada2) 
library(dplyr)
library(gridExtra)
library(SummarizedExperiment)
library(tidyr)
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
library(reshape2)
library(microbiomeutilities)
library(ggpubr)
library(ggpmisc)
library(nlme)
library(lmerTest)

##### Load important starting objects or load a workspace #####
#Start with the existing objects from the phyloseq data with all experimental timepoints. 
#This follow from objects generated in the "16Sanalyses_script.R". 
psexp3000 <- readRDS("R_objects/psexp3000_gg2.rds")
dfexp3000 <- read.csv("psexp3000_metadat.csv")

#Or, load an existing workspace such as 
#load("AlphaBetaChangeWorkspace_6October2025.RData")


##### Pull metadata #####

psexp3000 #(rarefied to 3000 reads, excludes 18RLP weird animal)
#View(dfexp3000) #metadata pulled out from the sample data of the phylsoeq object, plus richness and a couple other things

psexp3000no31 <- subset_samples(psexp3000, day_num != "31", )
psexp3000no31 <- subset_samples(psexp3000no31, sample_names(psexp3000no31) != "21_RP_D28")
psexp3000no31 <- prune_taxa(taxa_sums(psexp3000no31)>1, psexp3000no31)
psexp3000
psexp3000no31 #Ok, that should now work. 

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

#Factor out any variables you want sorted in a specific way. 
dfexp3000no31$microb <- factor(dfexp3000no31$microb, levels=c("LO", "HI"))
dfexp3000no31$MT <- paste(dfexp3000no31$microb, dfexp3000no31$treatment, sep="_")
dfexp3000no31$MT <- factor(dfexp3000no31$MT, levels=c("LO_worm" ,"LO_control", "HI_worm", "HI_control"))


#Simplify this metadata so it's easier to handle. Right now, there are lots of columns I don't really need. 
dfmetadat <- select(dfexp3000no31, -c('run', 'number', 'plate', 'well', 'ear', 'sample_type', 'samp_date', 'out_of_bubble',
                            'DOB', 'moved_to_SB', 'moved_from_bubble', 'location', 'breeding', 'sample_group',
                            'description', 'LanePosition', 'PlatePosition', 'Nanostring_sample', 'metabolomics_sample', 
                            'cage_day'))
ncol(dfmetadat)

#
#   ##### Now add qPCR data to the dataframe so that we can just get everything in at once. #####
newqpcr <- read.csv("~/Documents/1U4U_qPCR/map_1U4U_Full16S_metadata_ByOrigID_wBactQuant.csv") #This is to avoid remaking below, with additional columns.  

colnames(newqpcr)
trimdf <- newqpcr[,c(1, 26, 30, 31)]
#View(trimdf) #That looks good. 
#rename so that it merges correctly. 

trimdf <- trimdf %>% rename("row" = "ID", 
                            "qPCR_input_DNA" = "input_DNA")
dfmetadat <- merge(dfmetadat, trimdf, by = "row", all.x = TRUE)
dfmetadat$copies_16S_per_ng <- trunc(dfmetadat$copies_16S_per_ng)
#

##### Summarize the richness, etc in each group at differnt timepoints #####

dfD1 <- dfmetadat[dfmetadat$day_num==1,]
dfD10 <- dfmetadat[dfmetadat$day_num==10,]
dfD14 <- dfmetadat[dfmetadat$day_num==14,]
dfD28 <- dfmetadat[dfmetadat$day_num==28,]
dfD41 <- dfmetadat[dfmetadat$day_num==41,]

summary(dfD28$Observed[dfD28$microb == "LO" & dfD28$treatment=="worm"])
sd(dfD28$Observed[dfD28$microb == "LO"& dfD28$treatment=="worm"])
#pretreatment lo_worm: 28.9 [± 11.89]
#1 dpi lo_worm: 31.8 [± 12.12] 
#10 dpi lo_worm: 64.7 [± 13.07] 
#14 dpi lo_worm: 60.19 [± 9.965869]
#28 dpi lo_worm: 63.6 [8.64209]

#1 dpi lo_sham: 43.83 +- 23.93672
#10 dpi lo_sham: 52.2 [± 15.28] 

#1 dpi hi_worm: 170.4 [± 33.2] 
#10 dpi hi_worm: 162.8 {+- 40.92928}
#14 dpi hi_worm: 164.9 [+- 33.022]
#28 dpi hi_worm: 99.6 [+- 24.81532]

#################### RESILIENCE CALCULATIONS - PRETREATMENT TO 28 AND 41 DPI #####
#   ##### Calculate beta diversity dissimilarity matrices #####
ps.rel <- microbiome::transform(psexp3000no31, transform="compositional")

#Make sure your objects are wht you want to be working with. 
pspre.rel
ps28.rel
ps41.rel

#Merge the phyloseq objects to get the preterat vs. each timepoint. 
pspre.28.rel <- merge_phyloseq(pspre.rel, ps28.rel)
pspre.41.rel <- merge_phyloseq(pspre.rel, ps41.rel)
#      ##### Now grab GUniFrac distances #####
obj <- ps.rel
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Generate the matrices you care about.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
head(d50)
#Start cleaning up the data. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#Now will need to break up the strings so that we can identify the samples and make sure we're getting the comparison of interest. 
#library(tidyr)
GU50.m.update <- GU50.m %>%
  separate(Var1, c("cage1", "ear1", "day1"), sep = "_", remove=FALSE)
GU50.m.update$mouse1 <- paste(GU50.m.update$cage1, GU50.m.update$ear1, sep="_")
GU50.m.update <- GU50.m.update %>%
  separate(Var2, c("cage2", "ear2", "day2"), sep = "_", remove=FALSE)
GU50.m.update$mouse2 <- paste(GU50.m.update$cage2, GU50.m.update$ear2, sep="_")

#Now keep the samples only that have pretreatment as one of the days. 
GU50sub1 <- GU50.m.update[GU50.m.update$day1 == "pre" &
                            GU50.m.update$mouse1 == GU50.m.update$mouse2,]
GU50sub2 <- GU50.m.update[GU50.m.update$day2 == "pre" &
                            GU50.m.update$mouse1 == GU50.m.update$mouse2,]
GU50sub2rename <- GU50sub2 %>%  #Rename the columns so that all the "pre" samples are all day1 so that I can merge them. 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1, 
    day1=day2, 
    day2=day1)

dftmp <- rbind(GU50sub1, GU50sub2rename)
#
GU50resildf <- dftmp[,c("value", "Var2")]
GU50resildf <- GU50resildf %>% rename("row" = "Var2", 
                                      "GU50" = "value")
#

#   ##### Bind beta diversity with alpha diversity #####
#Main df is dfmetadat. 
colnames(dfmetadat)
head(GU50resildf)
#Merge the two betadiv dfs together. 
#
#Now merge with the main metadat. 
dfmetadatwhole <- merge(dfmetadat, GU50resildf, by = "row", all.x=TRUE)
#
#Make all NAs for pretreatment into 0. 
dfmetadatwhole$GU50[is.na(dfmetadatwhole$GU50)] <- 0
#

##### REARRANGE DF TO CALCULATE CHANGE in all metrics btwn pretreatment and all other timepoints. #####
dfpre <- subset(dfmetadatwhole, day=="D-5") #df of only D-5 samples, from df where I'd already added
# alpha div info
colnames(dfpre)
dfpre <- dfpre %>% rename("obs.pre" = "Observed",
                        "shan.pre" = "Shannon",
                        "faiths.pre" = "PD", 
                        "evenness.pre" = "evenness",
                        "copies_16S.pre" = "copies_16S",
                        "copies_16S_per_ng.pre" = "copies_16S_per_ng", 
                        "GU50.pre" = "GU50")


#Now pull out the other days and create subdf's for each day. Keep the mouse column for future merging.
dfD1 <- subset(dfmetadatwhole,day=="D1")[,c(6,18,20,21,23,25,27:29)]
dfD1 <- dfD1 %>% rename("obs.D1" = "Observed",
                        "shan.D1" = "Shannon",
                        "faiths.D1" = "PD", 
                        "evenness.D1" = "evenness",
                        "copies_16S.D1" = "copies_16S",
                        "copies_16S_per_ng.D1" = "copies_16S_per_ng", 
                        "GU50.D1" = "GU50")
                        
dfD10 <- subset(dfmetadatwhole,day=="D10")[,c(6,18,20,21,23,25,27:29)]
dfD10 <- dfD10 %>% rename("obs.D10" = "Observed",
                        "shan.D10" = "Shannon",
                        "faiths.D10" = "PD", 
                        "evenness.D10" = "evenness",
                        "copies_16S.D10" = "copies_16S",
                        "copies_16S_per_ng.D10" = "copies_16S_per_ng", 
                        "GU50.D10" = "GU50")

dfD14 <- subset(dfmetadatwhole,day=="D14")[,c(6,18,20,21,23,25,27:29)] #df of only D14 unique info.
dfD14 <- dfD14 %>% rename("obs.D14" = "Observed",
                          "shan.D14" = "Shannon",
                          "faiths.D14" = "PD", 
                          "evenness.D14" = "evenness",
                          "copies_16S.D14" = "copies_16S",
                          "copies_16S_per_ng.D14" = "copies_16S_per_ng", 
                          "GU50.D14" = "GU50")

dfD28 <- subset(dfmetadatwhole,day=="D28")[,c(6,18,20,21,23,25,27:29)]
dfD28 <- dfD28 %>% rename("obs.D28" = "Observed",
                          "shan.D28" = "Shannon",
                          "faiths.D28" = "PD", 
                          "evenness.D28" = "evenness",
                          "copies_16S.D28" = "copies_16S",
                          "copies_16S_per_ng.D28" = "copies_16S_per_ng", 
                          "GU50.D28" = "GU50")

dfD41 <- subset(dfmetadatwhole,day=="D41")[,c(6,18,20,21,23,25,27:29)]
dfD41 <- dfD41 %>% rename("obs.D41" = "Observed",
                          "shan.D41" = "Shannon",
                          "faiths.D41" = "PD", 
                          "evenness.D41" = "evenness",
                          "copies_16S.D41" = "copies_16S",
                          "copies_16S_per_ng.D41" = "copies_16S_per_ng", 
                          "GU50.D41" = "GU50")


dftmp <- merge(dfpre, dfD1, by="mouse", all=TRUE)
dftmp <- merge(dftmp, dfD10, by="mouse", all=TRUE)
dftmp <- merge(dftmp, dfD14, by="mouse", all=TRUE)
dftmp <- merge(dftmp, dfD28, by="mouse", all=TRUE)
dfmerge <- merge(dftmp, dfD41, by="mouse", all=TRUE)
#View(dfmerge) #Make sure it all came out ok. 
#

#   ##### Calculate change in all metrics between pretreatment and each other timepoint #####
#Add alpha div difference cols. 
dfmerge$start <- 0
dfmerge$obspre1 <- dfmerge$obs.D1-dfmerge$obs.pre
dfmerge$obspre10 <- dfmerge$obs.D10-dfmerge$obs.pre
dfmerge$obspre14 <- dfmerge$obs.D14-dfmerge$obs.pre
dfmerge$obspre28 <- dfmerge$obs.D28-dfmerge$obs.pre
dfmerge$obspre41 <- dfmerge$obs.D41-dfmerge$obs.pre

dfmerge$evennesspre1 <- dfmerge$evenness.D1-dfmerge$evenness.pre
dfmerge$evennesspre10 <- dfmerge$evenness.D10-dfmerge$evenness.pre
dfmerge$evennesspre14 <- dfmerge$evenness.D14-dfmerge$evenness.pre
dfmerge$evennesspre28 <- dfmerge$evenness.D28-dfmerge$evenness.pre
dfmerge$evennesspre41 <- dfmerge$evenness.D41-dfmerge$evenness.pre

dfmerge$shanpre1 <- dfmerge$shan.D1-dfmerge$shan.pre
dfmerge$shanpre10 <- dfmerge$shan.D10-dfmerge$shan.pre
dfmerge$shanpre14 <- dfmerge$shan.D14-dfmerge$shan.pre
dfmerge$shanpre28 <- dfmerge$shan.D28-dfmerge$shan.pre
dfmerge$shanpre41 <- dfmerge$shan.D41-dfmerge$shan.pre

dfmerge$faithspre1 <- dfmerge$faiths.D1-dfmerge$faiths.pre
dfmerge$faithspre10 <- dfmerge$faiths.D10-dfmerge$faiths.pre
dfmerge$faithspre14 <- dfmerge$faiths.D14-dfmerge$faiths.pre
dfmerge$faithspre28 <- dfmerge$faiths.D28-dfmerge$faiths.pre
dfmerge$faithspre41 <- dfmerge$faiths.D41-dfmerge$faiths.pre

dfmerge$GU50pre1 <- dfmerge$GU50.D1-dfmerge$GU50.pre
dfmerge$GU50pre10 <- dfmerge$GU50.D10-dfmerge$GU50.pre
dfmerge$GU50pre14 <- dfmerge$GU50.D14-dfmerge$GU50.pre
dfmerge$GU50pre28 <- dfmerge$GU50.D28-dfmerge$GU50.pre
dfmerge$GU50pre41 <- dfmerge$GU50.D41-dfmerge$GU50.pre

dfmerge$copies_16Spre1 <- dfmerge$copies_16S.D1-dfmerge$copies_16S.pre
dfmerge$copies_16Spre10 <- dfmerge$copies_16S.D10-dfmerge$copies_16S.pre
dfmerge$copies_16Spre14 <- dfmerge$copies_16S.D14-dfmerge$copies_16S.pre
dfmerge$copies_16Spre28 <- dfmerge$copies_16S.D28-dfmerge$copies_16S.pre
dfmerge$copies_16Spre41 <- dfmerge$copies_16S.D41-dfmerge$copies_16S.pre

dfmerge$copies_16S_per_ngpre1 <- dfmerge$copies_16S_per_ng.D1-dfmerge$copies_16S_per_ng.pre
dfmerge$copies_16S_per_ngpre10 <- dfmerge$copies_16S_per_ng.D10-dfmerge$copies_16S_per_ng.pre
dfmerge$copies_16S_per_ngpre14 <- dfmerge$copies_16S_per_ng.D14-dfmerge$copies_16S_per_ng.pre
dfmerge$copies_16S_per_ngpre28 <- dfmerge$copies_16S_per_ng.D28-dfmerge$copies_16S_per_ng.pre
dfmerge$copies_16S_per_ngpre41 <- dfmerge$copies_16S_per_ng.D41-dfmerge$copies_16S_per_ng.pre

#Add fractional change in ASV richness. 
dfmerge$obspre1perc <- dfmerge$obspre1/dfmerge$obs.pre
dfmerge$obspre10perc <- dfmerge$obspre10/dfmerge$obs.pre
dfmerge$obspre14perc <- dfmerge$obspre14/dfmerge$obs.pre
dfmerge$obspre28perc <- dfmerge$obspre28/dfmerge$obs.pre
dfmerge$obspre41perc <- dfmerge$obspre41/dfmerge$obs.pre

#   ##### Melt df's and put back together. #####
#Really though, we should melt the df to get it into a more graphable form over time. 
library(reshape2)
colnames(dfmerge)

#Melt Richness Change 
dfmeltobs <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "obspre1", "obspre10", "obspre14", "obspre28", "obspre41"), 
                  value.name="ChangeASVs", variable.name="time")
head(dfmeltobs)
#Melt Evenness Change
dfmeltevenness <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "evennesspre1", "evennesspre10", "evennesspre14", "evennesspre28", "evennesspre41"), 
                       value.name="ChangeEvenness", variable.name="time")
#Melt Shannon Change 
dfmeltshan <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "shanpre1", "shanpre10", "shanpre14", "shanpre28", "shanpre41"), 
                  value.name="ChangeShannon", variable.name="time")
#Melt Faith's Change 
dfmeltfaiths <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "faithspre1", "faithspre10", "faithspre14", "faithspre28", "faithspre41"), 
                  value.name="ChangeFaiths", variable.name="time")
#Melt Bray-Curtis Change 
dfmeltbray <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "braypre1", "braypre10", "braypre14", "braypre28", "braypre41"), 
                  value.name="ChangeBray", variable.name="time")
#Melt GUniFrac Change 
dfmeltGU50 <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "GU50pre1", "GU50pre10", "GU50pre14", "GU50pre28", "GU50pre41"), 
                  value.name="ChangeGUniFrac", variable.name="time")
#Melt copies 16S Change 
dfmeltcopies_16S <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "copies_16Spre1", "copies_16Spre10", "copies_16Spre14", "copies_16Spre28", "copies_16Spre41"), 
                  value.name="ChangeCopies16S", variable.name="time")
#Melt Richness Change 
dfmeltcopies_16S_per_ng <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "copies_16S_per_ngpre1", "copies_16S_per_ngpre10", "copies_16S_per_ngpre14", "copies_16S_per_ngpre28", "copies_16S_per_ngpre41"), 
                  value.name="ChangeCopies16SperNG", variable.name="time")
#Melt for fractional ASV richness change.
dfpercResil <- melt(dfmerge, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "obspre1perc", "obspre10perc", "obspre14perc", "obspre28perc", "obspre41perc"), 
                    value.name="ChangePercASVs", variable.name="time")

#Stick the change in all the metrics together to form 1 df. 
dfrawResil <- cbind(dfmeltobs, 
                    "ChangeEvenness" = dfmeltevenness$ChangeEvenness, 
                    "ChangeShannon" = dfmeltshan$ChangeShannon, 
                    "ChangeFaiths" = dfmeltfaiths$ChangeFaiths, 
                    "ChangeBray" = dfmeltbray$ChangeBray, 
                    "ChangeGUniFrac" = dfmeltGU50$ChangeGUniFrac, 
                    "ChangeCopies16S" = dfmeltcopies_16S$ChangeCopies16S, 
                    "ChangeCopies16SperNG" = dfmeltcopies_16S_per_ng$ChangeCopies16SperNG, 
                    "ChangePercASVs" = dfpercResil$ChangePercASVs)

#Create a new column that will eventually be numeric day
dfrawResil$day_num <- as.character(dfrawResil$time)
#Replace the values to make numeric
dfrawResil <- dfrawResil %>% 
  mutate(day_num = replace(day_num, day_num == 'start', -5)) %>%
  mutate(day_num = replace(day_num, day_num == 'obspre1', 1)) %>%
  mutate(day_num = replace(day_num, day_num == 'obspre10', 10)) %>%
  mutate(day_num = replace(day_num, day_num == 'obspre14', 14)) %>%
  mutate(day_num = replace(day_num, day_num == 'obspre28', 28)) %>%
  mutate(day_num = replace(day_num, day_num == 'obspre41', 41))
#Now make the column numeric. 
dfrawResil$day_num <- as.numeric(dfrawResil$day_num)

View(dfrawResil)
nrow(dfrawResil)
#
#Remove the rows for D28 for all HI. Lacking data.  
dfrawResil1 <- dfrawResil[!(dfrawResil$day_num == "28" & dfrawResil$MT == "HI_control"),]
dfrawResil2 <- dfrawResil1[!(dfrawResil1$day_num == "41" & dfrawResil1$mouse == "21_LP"),]
nrow(dfrawResil2)
#
#Ok, now drop all the rows that don't have any data for all metrics. 
dfResil <- dfrawResil2 %>% tidyr::drop_na(ChangeASVs)
nrow(dfResil) #final number of rows is 303
colnames(dfResil)
#

#   ##### Add MDT column for colors #####
dfResil$MDT <- paste(dfResil$microb, dfResil$day_num, dfResil$treatment, sep="_")
dfResil$MDT <- factor(dfResil$MDT, levels=c('LO_-5_worm', 'LO_1_worm', 'LO_10_worm', 'LO_14_worm', 
                                            'LO_28_worm', 'LO_31_worm', 'LO_41_worm', 'HI_-5_worm', 
                                            'HI_1_worm', 'HI_10_worm', 'HI_14_worm', 'HI_28_worm', 
                                            'HI_31_worm', 'HI_41_worm', 'LO_-5_control', 'LO_1_control', 
                                            'LO_10_control', 'LO_14_control', 'LO_28_control', 'LO_31_control', 
                                            'LO_41_control', 'HI_-5_control', 'HI_1_control', 'HI_10_control', 
                                            'HI_14_control', 'HI_28_control', 'HI_31_control', 'HI_41_control'))

#Visualize pattern of change over time. 


##### Example plots of change from pretreatment to every other timepoint #####
#   ##### Quick and dirty plot of each change column just to see what general trends look like ##### 
plota <- ggplot(dfResil, aes(x=day_num, y=ChangeASVs)) + geom_point(aes(color=MT)) +theme_bw()

plotb <- ggplot(dfResil, aes(x=day_num, y=ChangeEvenness)) + geom_point(aes(color=MT)) +theme_bw()
  
plotc <- ggplot(dfResil, aes(x=day_num, y=ChangeFaiths)) + geom_point(aes(color=MT)) +theme_bw()

grid.arrange(plota, plotb, plotc)

å#Real quick, what is going on with the HI_control individuals??
dftmp <- dfResil[dfResil$MT == "HI_control", ]
ggplot(dftmp, aes(x=day_num, y=ChangePercASVs)) + geom_point(aes(color=mouse)) +theme_bw()
#How about just LO_control? 
dftmp <- dfResil[dfResil$MT == "LO_control", ]
ggplot(dftmp, aes(x=day_num, y=ChangeASVs)) + geom_point(aes(color=mouse)) +theme_bw()
#Note that cage 12 changed a lot, and the other two cages were relatively unchanged. This is lost in the averages
# but important since it's one cage that's pulling the averages up. 

#How about looking at each of the groups on separates panes?
ggplot(dfResil, aes(x=day_num, y=ChangePercASVs)) + geom_point(aes(color=MT)) +
  theme_bw() +
  facet_wrap(~MT) +
  ggtitle("Percent Change Richness faceted by Microbiome/Treatment subgroups")

#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/ChangePercRichness_dot_FacetedbyMT.pdf", units="in", height=9, width=10)

#
#   ##### Example line plots over all timepints with error bars #####
#      ##### ASV Richness Create means and sd #####
dfResilhitreat <- dfResil[dfResil$microb=="HI" & dfResil$treatment=="worm",]
dfResilhitreat <- dfResilhitreat %>% tidyr::drop_na(ChangeASVs)
dfResilhictrl <- dfResil[dfResil$microb=="HI" & dfResil$treatment=="control",]
dfResilhictrl <- dfResilhictrl %>% tidyr::drop_na(ChangeASVs)
dfResillotreat <- dfResil[dfResil$microb=="LO" & dfResil$treatment=="worm",]
dfResillotreat <- dfResillotreat %>% tidyr::drop_na(ChangeASVs)
dfResilloctrl <- dfResil[dfResil$microb=="LO" & dfResil$treatment=="control",]
dfResilloctrl <- dfResilloctrl %>% tidyr::drop_na(ChangeASVs)


#Calculate the means of each group. 
dfResilhitreat.mean <- dfResilhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeASVs = mean(ChangeASVs, na.rm=TRUE))
dfResilhictrl.mean <- dfResilhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeASVs = mean(ChangeASVs, na.rm=TRUE))
dfResillotreat.mean <- dfResillotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeASVs = mean(ChangeASVs, na.rm=TRUE))
dfResilloctrl.mean <- dfResilloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeASVs = mean(ChangeASVs, na.rm=TRUE))

#Get the standard deviation
dfResilhitreat.sd <- dfResilhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeASVs, na.rm=TRUE))
dfResilhictrl.sd <- dfResilhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeASVs, na.rm=TRUE))
dfResillotreat.sd <- dfResillotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeASVs, na.rm=TRUE))
dfResilloctrl.sd <- dfResilloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeASVs, na.rm=TRUE))

#
dfResilhitreat.mean$SD <- dfResilhitreat.sd$SD
dfResilhictrl.mean$SD <- dfResilhictrl.sd$SD
dfResillotreat.mean$SD <- dfResillotreat.sd$SD
dfResilloctrl.mean$SD <- dfResilloctrl.sd$SD


# add sample size and calculate SEM 
hitreat.N <- as.data.frame(plyr::count(dfResilhitreat$day_num)) 
hictrl.N <- as.data.frame(plyr::count(dfResilhictrl$day_num)) 
lotreat.N <- as.data.frame(plyr::count(dfResillotreat$day_num)) 
loctrl.N <- as.data.frame(plyr::count(dfResilloctrl$day_num)) 

# add N to dataframes 
dfResilhitreat.mean$N <- hitreat.N$freq
dfResilhictrl.mean$N <- hictrl.N$freq
dfResillotreat.mean$N <- lotreat.N$freq
dfResilloctrl.mean$N <- loctrl.N$freq

# add SE of each group to the mean. 

dfResilhitreat.mean$SEH <- dfResilhitreat.mean$meanChangeASVs + (dfResilhitreat.mean$SD/(sqrt(dfResilhitreat.mean$N)))
dfResilhictrl.mean$SEH <- dfResilhictrl.mean$meanChangeASVs + (dfResilhictrl.mean$SD/(sqrt(dfResilhictrl.mean$N)))
dfResillotreat.mean$SEH <- dfResillotreat.mean$meanChangeASVs + (dfResillotreat.mean$SD/(sqrt(dfResillotreat.mean$N)))
dfResilloctrl.mean$SEH <- dfResilloctrl.mean$meanChangeASVs + (dfResilloctrl.mean$SD/(sqrt(dfResilloctrl.mean$N)))

dfResilhitreat.mean$SEL <- dfResilhitreat.mean$meanChangeASVs - (dfResilhitreat.mean$SD/(sqrt(dfResilhitreat.mean$N)))
dfResilhictrl.mean$SEL <- dfResilhictrl.mean$meanChangeASVs - (dfResilhictrl.mean$SD/(sqrt(dfResilhictrl.mean$N)))
dfResillotreat.mean$SEL <- dfResillotreat.mean$meanChangeASVs - (dfResillotreat.mean$SD/(sqrt(dfResillotreat.mean$N)))
dfResilloctrl.mean$SEL <- dfResilloctrl.mean$meanChangeASVs - (dfResilloctrl.mean$SD/(sqrt(dfResilloctrl.mean$N)))


#save these as ASVs-specific dfs. 
dfResilhitreat.mean.ASVs <- dfResilhitreat.mean
dfResilhictrl.mean.ASVs <- dfResilhictrl.mean 
dfResillotreat.mean.ASVs <- dfResillotreat.mean
dfResilloctrl.mean.ASVs <- dfResilloctrl.mean 

#         ##### Plot as average over time with error bars #####
MeanASVChangeplot <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfResilhitreat.mean.ASVs$meanChangeASVs, x=dfResilhitreat.mean.ASVs$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfResilhictrl.mean.ASVs$meanChangeASVs, x=dfResilhictrl.mean.ASVs$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfResilhitreat.mean.ASVs$SEL, ymax=dfResilhitreat.mean.ASVs$SEH, y=dfResilhitreat.mean.ASVs$meanChangeASVs, x=dfResilhitreat.mean.ASVs$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfResilhictrl.mean.ASVs$SEL, ymax=dfResilhictrl.mean.ASVs$SEH, y=dfResilhictrl.mean.ASVs$meanChangeASVs, x=dfResilhictrl.mean.ASVs$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfResilhitreat.mean.ASVs$meanChangeASVs, x=dfResilhitreat.mean.ASVs$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfResilhictrl.mean.ASVs$meanChangeASVs, x=dfResilhictrl.mean.ASVs$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dfResillotreat.mean.ASVs$meanChangeASVs, x=dfResillotreat.mean.ASVs$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfResilloctrl.mean.ASVs$meanChangeASVs, x=dfResilloctrl.mean.ASVs$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dfResillotreat.mean.ASVs$SEL, ymax=dfResillotreat.mean.ASVs$SEH, y=dfResillotreat.mean.ASVs$meanChangeASVs, x=dfResillotreat.mean.ASVs$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfResilloctrl.mean.ASVs$SEL, ymax=dfResilloctrl.mean.ASVs$SEH, y=dfResilloctrl.mean.ASVs$meanChangeASVs, x=dfResilloctrl.mean.ASVs$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dfResillotreat.mean.ASVs$meanChangeASVs, x=dfResillotreat.mean.ASVs$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfResilloctrl.mean.ASVs$meanChangeASVs, x=dfResilloctrl.mean.ASVs$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("Change in ASV richness since pretreatment") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average change in ASVs from pretreatment to each timepoint") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

MeanASVChangeplot
#
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/psexp3000_AverageChangeASVs_Resilience_lineplot_woD31.pdf", MeanASVChangeplot, units="in", height=5, width=6)
#


#      ##### GUniFrac Create means and sd #####
dfResilhitreat <- dfResil[dfResil$microb=="HI" & dfResil$treatment=="worm",]
dfResilhitreat <- dfResilhitreat %>% tidyr::drop_na(ChangeGUniFrac)
dfResilhictrl <- dfResil[dfResil$microb=="HI" & dfResil$treatment=="control",]
dfResilhictrl <- dfResilhictrl %>% tidyr::drop_na(ChangeGUniFrac)
dfResillotreat <- dfResil[dfResil$microb=="LO" & dfResil$treatment=="worm",]
dfResillotreat <- dfResillotreat %>% tidyr::drop_na(ChangeGUniFrac)
dfResilloctrl <- dfResil[dfResil$microb=="LO" & dfResil$treatment=="control",]
dfResilloctrl <- dfResilloctrl %>% tidyr::drop_na(ChangeGUniFrac)


#Calculate the means of each group. 
dfResilhitreat.mean <- dfResilhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeGUniFrac = mean(ChangeGUniFrac, na.rm=TRUE))
dfResilhictrl.mean <- dfResilhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeGUniFrac = mean(ChangeGUniFrac, na.rm=TRUE))
dfResillotreat.mean <- dfResillotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeGUniFrac = mean(ChangeGUniFrac, na.rm=TRUE))
dfResilloctrl.mean <- dfResilloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanChangeGUniFrac = mean(ChangeGUniFrac, na.rm=TRUE))

#Get the standard deviation
dfResilhitreat.sd <- dfResilhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeGUniFrac, na.rm=TRUE))
dfResilhictrl.sd <- dfResilhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeGUniFrac, na.rm=TRUE))
dfResillotreat.sd <- dfResillotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeGUniFrac, na.rm=TRUE))
dfResilloctrl.sd <- dfResilloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(ChangeGUniFrac, na.rm=TRUE))

#
dfResilhitreat.mean$SD <- dfResilhitreat.sd$SD
dfResilhictrl.mean$SD <- dfResilhictrl.sd$SD
dfResillotreat.mean$SD <- dfResillotreat.sd$SD
dfResilloctrl.mean$SD <- dfResilloctrl.sd$SD


# add sample size and calculate SEM 
hitreat.N <- as.data.frame(plyr::count(dfResilhitreat$day_num)) 
hictrl.N <- as.data.frame(plyr::count(dfResilhictrl$day_num)) 
lotreat.N <- as.data.frame(plyr::count(dfResillotreat$day_num)) 
loctrl.N <- as.data.frame(plyr::count(dfResilloctrl$day_num)) 

# add N to dataframes 
dfResilhitreat.mean$N <- hitreat.N$freq
dfResilhictrl.mean$N <- hictrl.N$freq
dfResillotreat.mean$N <- lotreat.N$freq
dfResilloctrl.mean$N <- loctrl.N$freq

# add SE of each group to the mean. 

dfResilhitreat.mean$SEH <- dfResilhitreat.mean$meanChangeGUniFrac + (dfResilhitreat.mean$SD/(sqrt(dfResilhitreat.mean$N)))
dfResilhictrl.mean$SEH <- dfResilhictrl.mean$meanChangeGUniFrac + (dfResilhictrl.mean$SD/(sqrt(dfResilhictrl.mean$N)))
dfResillotreat.mean$SEH <- dfResillotreat.mean$meanChangeGUniFrac + (dfResillotreat.mean$SD/(sqrt(dfResillotreat.mean$N)))
dfResilloctrl.mean$SEH <- dfResilloctrl.mean$meanChangeGUniFrac + (dfResilloctrl.mean$SD/(sqrt(dfResilloctrl.mean$N)))

dfResilhitreat.mean$SEL <- dfResilhitreat.mean$meanChangeGUniFrac - (dfResilhitreat.mean$SD/(sqrt(dfResilhitreat.mean$N)))
dfResilhictrl.mean$SEL <- dfResilhictrl.mean$meanChangeGUniFrac - (dfResilhictrl.mean$SD/(sqrt(dfResilhictrl.mean$N)))
dfResillotreat.mean$SEL <- dfResillotreat.mean$meanChangeGUniFrac - (dfResillotreat.mean$SD/(sqrt(dfResillotreat.mean$N)))
dfResilloctrl.mean$SEL <- dfResilloctrl.mean$meanChangeGUniFrac - (dfResilloctrl.mean$SD/(sqrt(dfResilloctrl.mean$N)))


#save these as GUniFrac-specific dfs. 
dfResilhitreat.mean.GUniFrac <- dfResilhitreat.mean
dfResilhictrl.mean.GUniFrac <- dfResilhictrl.mean 
dfResillotreat.mean.GUniFrac <- dfResillotreat.mean
dfResilloctrl.mean.GUniFrac <- dfResilloctrl.mean 

#         ##### Plot as average over time with error bars #####
MeanGUniFracChangeplot <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfResilhitreat.mean.GUniFrac$meanChangeGUniFrac, x=dfResilhitreat.mean.GUniFrac$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfResilhictrl.mean.GUniFrac$meanChangeGUniFrac, x=dfResilhictrl.mean.GUniFrac$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfResilhitreat.mean.GUniFrac$SEL, ymax=dfResilhitreat.mean.GUniFrac$SEH, y=dfResilhitreat.mean.GUniFrac$meanChangeGUniFrac, x=dfResilhitreat.mean.GUniFrac$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfResilhictrl.mean.GUniFrac$SEL, ymax=dfResilhictrl.mean.GUniFrac$SEH, y=dfResilhictrl.mean.GUniFrac$meanChangeGUniFrac, x=dfResilhictrl.mean.GUniFrac$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfResilhitreat.mean.GUniFrac$meanChangeGUniFrac, x=dfResilhitreat.mean.GUniFrac$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfResilhictrl.mean.GUniFrac$meanChangeGUniFrac, x=dfResilhictrl.mean.GUniFrac$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dfResillotreat.mean.GUniFrac$meanChangeGUniFrac, x=dfResillotreat.mean.GUniFrac$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfResilloctrl.mean.GUniFrac$meanChangeGUniFrac, x=dfResilloctrl.mean.GUniFrac$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dfResillotreat.mean.GUniFrac$SEL, ymax=dfResillotreat.mean.GUniFrac$SEH, y=dfResillotreat.mean.GUniFrac$meanChangeGUniFrac, x=dfResillotreat.mean.GUniFrac$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfResilloctrl.mean.GUniFrac$SEL, ymax=dfResilloctrl.mean.GUniFrac$SEH, y=dfResilloctrl.mean.GUniFrac$meanChangeGUniFrac, x=dfResilloctrl.mean.GUniFrac$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dfResillotreat.mean.GUniFrac$meanChangeGUniFrac, x=dfResillotreat.mean.GUniFrac$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfResilloctrl.mean.GUniFrac$meanChangeGUniFrac, x=dfResilloctrl.mean.GUniFrac$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("Change in GUniFrac distance since pretreatment") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average change in GUniFrac from pretreatment to each timepoint") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

MeanGUniFracChangeplot

ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/psexp3000_AverageChangeGUniFrac_Resilience_lineplot_woD31.pdf", MeanGUniFracChangeplot, units="in", height=5, width=6)
#


#   ##### Example boxplots for intervals of interest for Resilience #####
#      ##### ASV Richness #####
#Boxplot change in richness pretreatment to D1
pRichnesspre1 <- ggplot(dfD1resil, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD1resil$ChangeASVs), max(dfD1resil$ChangeASVs)+0.1*max(dfD1resil$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pretreatment-D1") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichnesspre1
#
#Plotting change in richness pretreatment to D10
pRichnesspre10 <- ggplot(dfD10resil, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD10resil$ChangeASVs), max(dfD10resil$ChangeASVs)+0.1*max(dfD10resil$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pre-D10") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichnesspre10
#
#Plotting change in richness pretreatment to D14
pRichnesspre14 <- ggplot(dfD14resil, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD14resil$ChangeASVs), max(dfD14resil$ChangeASVs)+0.1*max(dfD14resil$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pre-D14") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichnesspre14
#
#Boxplot change in richness pretreatment to D28
pRichnesspre28  <- ggplot(dfD28resil, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964")) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD28resil$ChangeASVs), max(dfD28resil$ChangeASVs)+0.1*max(dfD28resil$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pre-D28") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichnesspre28
#
#Boxplot change in richness pretreatment to D41
pRichnesspre41  <- ggplot(dfD41resil, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD41resil$ChangeASVs), max(dfD41resil$ChangeASVs)+0.1*max(dfD41resil$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pre-D41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichnesspre41
#

#         ##### Save the plots for Richness #####
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resilience/psexp3000_ChangeRichnesspre1_boxplot_woD31.pdf", pRichnesspre1, height = 6, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resilience/psexp3000_ChangeRichnesspre10_boxplot_woD31.pdf", pRichnesspre10, height = 6, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resilience/psexp3000_ChangeRichnesspre14_boxplot_woD31.pdf", pRichnesspre14, height = 6, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resilience/psexp3000_ChangeRichnesspre28_boxplot_woD31.pdf", pRichnesspre28, height = 6, width = 3, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resilience/psexp3000_ChangeRichnesspre41_boxplot_woD31.pdf", pRichnesspre41, height = 6, width = 3, units = "in")

#
#      ##### GUniFrac #####
#Boxplot change in richness pretreatment to D1
pGUniFracpre1 <- ggplot(dfD1resil, aes(x=MT, y=ChangeGUniFrac)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD1resil$ChangeGUniFrac), max(dfD1resil$ChangeGUniFrac)+0.1*max(dfD1resil$ChangeGUniFrac)) +
  ylab("Change in GUniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac distance pretreatment-D1") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFracpre1
#
#Plotting change in richness pre-D10
pGUniFracpre10 <- ggplot(dfD10resil, aes(x=MT, y=ChangeGUniFrac)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD10resil$ChangeGUniFrac), max(dfD10resil$ChangeGUniFrac)+0.1*max(dfD10resil$ChangeGUniFrac)) +
  ylab("Change in GUniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac distance pre-D10") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFracpre10
#
#Plotting change in richness pre-D14
pGUniFracpre14 <- ggplot(dfD14resil, aes(x=MT, y=ChangeGUniFrac)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD14resil$ChangeGUniFrac), max(dfD14resil$ChangeGUniFrac)+0.1*max(dfD14resil$ChangeGUniFrac)) +
  ylab("Change in GUniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac distance pre-D14") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFracpre14
#
#Boxplot change in richness pre-D28
pGUniFracpre28  <- ggplot(dfD28resil, aes(x=MT, y=ChangeGUniFrac)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964")) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD28resil$ChangeGUniFrac), max(dfD28resil$ChangeGUniFrac)+0.1*max(dfD28resil$ChangeGUniFrac)) +
  ylab("Change in GUniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac distance pre-D28") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFracpre28
#
#Boxplot change in richness pre-D41
pGUniFracpre41  <- ggplot(dfD41resil, aes(x=MT, y=ChangeGUniFrac)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD41resil$ChangeGUniFrac), max(dfD41resil$ChangeGUniFrac)+0.1*max(dfD41resil$ChangeGUniFrac)) +
  ylab("Change in GUniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac distance pre-D41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFracpre41
#

#         ##### Save the plots for GUniFrac distance #####
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resilience/psexp3000_GUniFracpre1_boxplot_woD31.pdf", pGUniFracpre1, height = 6, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resilience/psexp3000_GUniFracpre10_boxplot_woD31.pdf", pGUniFracpre10, height = 6, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resilience/psexp3000_GUniFracpre14_boxplot_woD31.pdf", pGUniFracpre14, height = 6, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resilience/psexp3000_GUniFracpre28_boxplot_woD31.pdf", pGUniFracpre28, height = 6, width = 3, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resilience/psexp3000_GUniFracpre41_boxplot_woD31.pdf", pGUniFracpre41, height = 6, width = 3, units = "in")
#
#



##### Example stats for pretreatment to 28 and 41 dpi #####
#   ##### Break Resilience dataset up into timepoints #####
dfD1resil <- dfResil[dfResil$day_num==1,]
dfD10resil <- dfResil[dfResil$day_num==10,]
dfD14resil <- dfResil[dfResil$day_num==14,]
dfD28resil <- dfResil[dfResil$day_num==28,]
dfD41resil <- dfResil[dfResil$day_num==41,]
#
#      ##### Run wilcoxons on each timepoint of interest for change in richness.  #####

#Testing change pretreatment to D28. 
df1 <- dfD28resil[dfD28resil$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 75, p-value = 0.001224 
df2 <- dfD28resil[dfD28resil$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 7.5, p-value = 0.003924
#Can't compute HI vs. HI bc no control samples. 
#
ps <- c(0.001224 , 0.003924)
p.adjust(ps, "holm")
## 0.002448, 0.003924


#Testing change pretreatment to D41. 
df1 <- dfD41resil[dfD41resil$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 161, p-value = 0.001973 
df2 <- dfD41resil[dfD41resil$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 22, p-value = 0.01217
df3 <- dfD41resil[dfD41resil$microb == "HI",]
wilcox.test(ChangeASVs ~ treatment, data=df3)
#W = 8.5, p-value = 0.5394
#
ps <- c(0.001973, 0.01217, 0.5394)
p.adjust(ps, "holm")
## 0.005919, 0.024340, 0.539400


#


#      ##### Run wilcoxons on each timepoint of interest for change in GU50.  #####

#df1 is always worm-treated vs. worm-treated. df2 is Lo vs. Lo. df3 is hi vs. hi
df1 <- dfD1resil[dfD1resil$treatment == "worm",]
wilcox.test(ChangeGUniFrac ~ microb, data=df1) 
#W = 148, p-value = 0.00408
df2 <- dfD1resil[dfD1resil$microb == "LO",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df2)
#W = 96, p-value = 1
df3 <- dfD1resil[dfD1resil$microb == "HI",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df3)
#W = 41, p-value = 0.712
#
ps <- c(0.00408, 1, 0.712)
p.adjust(ps, "holm")
##0.01224 1.00000 1.00000


#Testing change pretreatment to D10. 
df1 <- dfD10resil[dfD10resil$treatment == "worm",]
wilcox.test(ChangeGUniFrac ~ microb, data=df1)
#W = 314, p-value = 0.6096
df2 <- dfD10resil[dfD10resil$microb == "LO",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df2)
#W = 40, p-value = 0.02343
df3 <- dfD10resil[dfD10resil$microb == "HI",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df3)
#W = 52, p-value = 0.1952
#
ps <- c(0.6096, 0.02343, 0.1952)
p.adjust(ps, "holm")
##0.60960 0.07029 0.39040
#

#Testing change pretreatment to D14. 
df1 <- dfD14resil[dfD14resil$treatment == "worm",]
wilcox.test(ChangeGUniFrac ~ microb, data=df1)
#W = 404, p-value = 8.398e-05
df2 <- dfD14resil[dfD14resil$microb == "LO",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df2)
#W = 25, p-value = 0.002777
df3 <- dfD14resil[dfD14resil$microb == "HI",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df3)
#W = 50, p-value = 0.0485
#
ps <- c(0.00008398, 0.002777, 0.0485)
p.adjust(ps, "holm")
## 0.00025194 0.00555400 0.04850000


#Testing change pretreatment to D28. 
df1 <- dfD28resil[dfD28resil$treatment == "worm",]
wilcox.test(ChangeGUniFrac ~ microb, data=df1)
#W = 64, p-value = 0.01935
df2 <- dfD28resil[dfD28resil$microb == "LO",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df2)
#W = 13, p-value = 0.01098
#Can't compute HI vs. HI bc no control samples. 
#
ps <- c(0.01935 , 0.01098)
p.adjust(ps, "holm")
## 0.02196 0.02196


#Testing change pretreatment to D41. 
df1 <- dfD41resil[dfD41resil$treatment == "worm",]
wilcox.test(ChangeGUniFrac ~ microb, data=df1)
#W = 155, p-value = 0.003225
df2 <- dfD41resil[dfD41resil$microb == "LO",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df2)
#W = 54, p-value = 0.4451
df3 <- dfD41resil[dfD41resil$microb == "HI",]
wilcox.test(ChangeGUniFrac ~ treatment, data=df3)
#W = 12, p-value = 1
#
ps <- c(0.003225, 0.4451, 1)
p.adjust(ps, "holm")
## 0.009675 0.890200 1.000000
#


##################### RESISTANCE CALCULATIONS - CHANGE BETWEEN ADJACENT TIMEPOINTS. #####


#   ##### Calculate change in alpha div metrics and qPCR data between adjacent timepoints -- Resistance. #####
#Duplicate the merged DF to make resistance calculations (adjacent timepoints)
dfresistance <- dfmerge[,1:70]

#Add alpha div difference cols. 
dfresistance$obs1 <- dfresistance$obs.D1-dfresistance$obs.pre
dfresistance$obs10 <- dfresistance$obs.D10-dfresistance$obs.D1
dfresistance$obs14 <- dfresistance$obs.D14-dfresistance$obs.D10
dfresistance$obs28 <- dfresistance$obs.D28-dfresistance$obs.D14
dfresistance$obs41 <- dfresistance$obs.D41-dfresistance$obs.D28
dfresistance$obs1441 <- dfresistance$obs.D41-dfresistance$obs.D14
dfresistance$obs1428 <- dfresistance$obs.D28-dfresistance$obs.D14
dfresistance$obs2841 <- dfresistance$obs.D41-dfresistance$obs.D28

#Now for evenness
dfresistance$evenness1 <- dfresistance$evenness.D1-dfresistance$evenness.pre
dfresistance$evenness10 <- dfresistance$evenness.D10-dfresistance$evenness.D1
dfresistance$evenness14 <- dfresistance$evenness.D14-dfresistance$evenness.D10
dfresistance$evenness28 <- dfresistance$evenness.D28-dfresistance$evenness.D14
dfresistance$evenness41 <- dfresistance$evenness.D41-dfresistance$evenness.D28
dfresistance$evenness1441 <- dfresistance$evenness.D41-dfresistance$evenness.D14
dfresistance$evenness1428 <- dfresistance$evenness.D28-dfresistance$evenness.D14
dfresistance$evenness2841 <- dfresistance$evenness.D41-dfresistance$evenness.D28

#Add alpha div difference cols. 
dfresistance$shan1 <- dfresistance$shan.D1-dfresistance$shan.pre
dfresistance$shan10 <- dfresistance$shan.D10-dfresistance$shan.D1
dfresistance$shan14 <- dfresistance$shan.D14-dfresistance$shan.D10
dfresistance$shan28 <- dfresistance$shan.D28-dfresistance$shan.D14
dfresistance$shan41 <- dfresistance$shan.D41-dfresistance$shan.D28
dfresistance$shan1441 <- dfresistance$shan.D41-dfresistance$shan.D14
dfresistance$shan1428 <- dfresistance$shan.D28-dfresistance$shan.D14
dfresistance$shan2841 <- dfresistance$shan.D41-dfresistance$shan.D28

#Add alpha div difference cols. 
dfresistance$faiths1 <- dfresistance$faiths.D1-dfresistance$faiths.pre
dfresistance$faiths10 <- dfresistance$faiths.D10-dfresistance$faiths.D1
dfresistance$faiths14 <- dfresistance$faiths.D14-dfresistance$faiths.D10
dfresistance$faiths28 <- dfresistance$faiths.D28-dfresistance$faiths.D14
dfresistance$faiths41 <- dfresistance$faiths.D41-dfresistance$faiths.D28
dfresistance$faiths1441 <- dfresistance$faiths.D41-dfresistance$faiths.D14
dfresistance$faiths1428 <- dfresistance$faiths.D28-dfresistance$faiths.D14
dfresistance$faiths2841 <- dfresistance$faiths.D41-dfresistance$faiths.D28

#Add alpha div difference cols. 
dfresistance$copies_16S1 <- dfresistance$copies_16S.D1-dfresistance$copies_16S.pre
dfresistance$copies_16S10 <- dfresistance$copies_16S.D10-dfresistance$copies_16S.D1
dfresistance$copies_16S14 <- dfresistance$copies_16S.D14-dfresistance$copies_16S.D10
dfresistance$copies_16S28 <- dfresistance$copies_16S.D28-dfresistance$copies_16S.D14
dfresistance$copies_16S41 <- dfresistance$copies_16S.D41-dfresistance$copies_16S.D28
dfresistance$copies_16S1441 <- dfresistance$copies_16S.D41-dfresistance$copies_16S.D14
dfresistance$copies_16S1428 <- dfresistance$copies_16S.D28-dfresistance$copies_16S.D14
dfresistance$copies_16S2841 <- dfresistance$copies_16S.D41-dfresistance$copies_16S.D28

#Add alpha div difference cols. 
dfresistance$copies_16S_per_ng1 <- dfresistance$copies_16S_per_ng.D1-dfresistance$copies_16S_per_ng.pre
dfresistance$copies_16S_per_ng10 <- dfresistance$copies_16S_per_ng.D10-dfresistance$copies_16S_per_ng.D1
dfresistance$copies_16S_per_ng14 <- dfresistance$copies_16S_per_ng.D14-dfresistance$copies_16S_per_ng.D10
dfresistance$copies_16S_per_ng28 <- dfresistance$copies_16S_per_ng.D28-dfresistance$copies_16S_per_ng.D14
dfresistance$copies_16S_per_ng41 <- dfresistance$copies_16S_per_ng.D41-dfresistance$copies_16S_per_ng.D28
dfresistance$copies_16S_per_ng1441 <- dfresistance$copies_16S_per_ng.D41-dfresistance$copies_16S_per_ng.D14
dfresistance$copies_16S_per_ng1428 <- dfresistance$copies_16S_per_ng.D28-dfresistance$copies_16S_per_ng.D14
dfresistance$copies_16S_per_ng2841 <- dfresistance$copies_16S_per_ng.D41-dfresistance$copies_16S_per_ng.D28

#      ##### Create the final df for resistance starting at pretreatment before adding beta diversity. ######

#Melt the dfs so that we can then stack columns together. 
dfmeltobs <- melt(dfresistance, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "obs1", "obs10", "obs14", "obs28", "obs41"), 
                  value.name="ChangeASVs", variable.name="time")
head(dfmeltobs)
#Melt Evenness Change
dfmeltevenness <- melt(dfresistance, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "evenness1", "evenness10", "evenness14", "evenness28", "evenness41"), 
                       value.name="ChangeEvenness", variable.name="time")
#Melt Shannon Change 
dfmeltshan <- melt(dfresistance, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "shan1", "shan10", "shan14", "shan28", "shan41"), 
                   value.name="ChangeShannon", variable.name="time")
#Melt Faith's Change 
dfmeltfaiths <- melt(dfresistance, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "faiths1", "faiths10", "faiths14", "faiths28", "faiths41"), 
                     value.name="ChangeFaiths", variable.name="time")
#Melt copies 16S Change 
dfmeltcopies_16S <- melt(dfresistance, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "copies_16S1", "copies_16S10", "copies_16S14", "copies_16S28", "copies_16S41"), 
                         value.name="ChangeCopies16S", variable.name="time")
#Melt Richness Change 
dfmeltcopies_16S_per_ng <- melt(dfresistance, id.vars=c(1,2,5,7:10,17,24), measure.vars=c("start", "copies_16S_per_ng1", "copies_16S_per_ng10", "copies_16S_per_ng14", "copies_16S_per_ng28", "copies_16S_per_ng41"), 
                                value.name="ChangeCopies16SperNG", variable.name="time")

#Stick the change in all the metrics together to form 1 df. 
dfrawResist <- cbind(dfmeltobs, 
                    "ChangeEvenness" = dfmeltevenness$ChangeEvenness, 
                    "ChangeShannon" = dfmeltshan$ChangeShannon, 
                    "ChangeFaiths" = dfmeltfaiths$ChangeFaiths, 
                    "ChangeCopies16S" = dfmeltcopies_16S$ChangeCopies16S, 
                    "ChangeCopies16SperNG" = dfmeltcopies_16S_per_ng$ChangeCopies16SperNG)

#Create a new column that will eventually be numeric day
dfrawResist$day_num <- as.character(dfrawResist$time)
#Replace the values to make numeric
dfrawResist <- dfrawResist %>% 
  mutate(day_num = replace(day_num, day_num == 'start', -5)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs1', 1)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs10', 10)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs14', 14)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs28', 28)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs41', 41))
#Now make the column numeric. 
dfrawResist$day_num <- as.numeric(dfrawResist$day_num)

View(dfrawResist)
nrow(dfrawResist)
#
#Remove the rows for D28 for all HI. Lacking data.  
dfrawResist1 <- dfrawResist[!(dfrawResist$day_num == "28" & dfrawResist$MT == "HI_control"),]
dfrawResist2 <- dfrawResist1[!(dfrawResist1$day_num == "41" & dfrawResist1$mouse == "21_LP"),]
nrow(dfrawResist2)
#
#Ok, now drop all the rows that don't have any data for all metrics. 
dfResist <- dfrawResist2 %>% tidyr::drop_na(ChangeASVs)
nrow(dfResist) #final number of rows is 288 for Resistance. More dropped than resilience bc more depends on D28.
colnames(dfResist)
#
# Drop the existing "row" column and create a new column that will make it possible to merge with beta diversity metrics. 
dfResist$day <- paste0("D", dfResist$day_num)
dfResist <- dfResist %>% mutate(day = replace(day, day == 'D-5', 'pre'))
dfResist$row <- paste(dfResist$mouse, dfResist$day, sep = "_")

#
#      ##### Now calculate beta diversity distances for adjacent timepoints. ##### 
ps.rel <- microbiome::transform(psexp3000no31, transform="compositional")

#Make sure your objects are wht you want to be working with. 
pspre.rel
ps1.rel
ps10.rel
ps14.rel
ps28.rel
ps41.rel

#Merge the phyloseq objects to get the preterat vs. each timepoint. 
pspre.1.rel <- merge_phyloseq(pspre.rel, ps1.rel)
ps1.10.rel <- merge_phyloseq(ps1.rel, ps10.rel)
ps10.14.rel <- merge_phyloseq(ps10.rel, ps14.rel)
ps14.28.rel <- merge_phyloseq(ps14.rel, ps28.rel)
ps28.41.rel <- merge_phyloseq(ps28.rel, ps41.rel)

#         ##### Now do GUniFrac. #####
#            ##### Pretreatment to D1 #####
pspre.1.rel <- merge_phyloseq(pspre.rel, ps1.rel)

obj <- pspre.1.rel
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Generate the matrices you care about.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
head(d50)
#Start cleaning up the data. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#Now will need to break up the strings so that we can identify the samples and make sure we're getting the comparison of interest. 
#library(tidyr)
GU50.m.update <- GU50.m %>%
  separate(Var1, c("cage1", "ear1", "day1"), sep = "_", remove=FALSE)
GU50.m.update$mouse1 <- paste(GU50.m.update$cage1, GU50.m.update$ear1, sep="_")
GU50.m.update <- GU50.m.update %>%
  separate(Var2, c("cage2", "ear2", "day2"), sep = "_", remove=FALSE)
GU50.m.update$mouse2 <- paste(GU50.m.update$cage2, GU50.m.update$ear2, sep="_")

GU50.sd.1 <- GU50.sd[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                       GU50.m.update$day1=="pre" &
                       GU50.m.update$day2=="D1",]
GU50.sd.1 <- GU50.sd.1 %>% drop_na(value)
#
GU50.sd.2 <- GU50.m.update[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                             GU50.m.update$day1=="D1" &
                             GU50.m.update$day2=="pre",]
GU50.sd.2 <- GU50.sd.2 %>% drop_na(value)
#
GU50.sd.2rename <- GU50.sd.2 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1, 
    day1=day2, 
    day2=day1)
GU50.sd.pre1 <- rbind(GU50.sd.1, GU50.sd.2rename)
GU50.sd.pre1 <- GU50.sd.pre1 %>%  #Rename the rownames so that you can merge with metadata df.  
  dplyr::rename(row=Var2)


#
#            ##### from D1 to D10 #####
ps1.10.rel <- merge_phyloseq(ps1.rel, ps10.rel)

obj <- ps1.10.rel
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Generate the matrices you care about.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
head(d50)
#Start cleaning up the data. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#Now will need to break up the strings so that we can identify the samples and make sure we're getting the comparison of interest. 
#library(tidyr)
GU50.m.update <- GU50.m %>%
  separate(Var1, c("cage1", "ear1", "day1"), sep = "_", remove=FALSE)
GU50.m.update$mouse1 <- paste(GU50.m.update$cage1, GU50.m.update$ear1, sep="_")
GU50.m.update <- GU50.m.update %>%
  separate(Var2, c("cage2", "ear2", "day2"), sep = "_", remove=FALSE)
GU50.m.update$mouse2 <- paste(GU50.m.update$cage2, GU50.m.update$ear2, sep="_")

GU50.sd.1 <- GU50.sd[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                       GU50.m.update$day1=="D1" &
                       GU50.m.update$day2=="D10",]
GU50.sd.1 <- GU50.sd.1 %>% drop_na(value)
#
GU50.sd.2 <- GU50.m.update[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                             GU50.m.update$day1=="D10" &
                             GU50.m.update$day2=="D1",]
GU50.sd.2 <- GU50.sd.2 %>% drop_na(value)
#
GU50.sd.2rename <- GU50.sd.2 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1, 
    day1=day2, 
    day2=day1)
GU50.sd.110 <- rbind(GU50.sd.1, GU50.sd.2rename)
GU50.sd.110 <- GU50.sd.110 %>%  #Rename the rownames so that you can merge with metadata df.  
  dplyr::rename(row=Var2)

#
#
#            ##### from D10 to D14 #####

ps10.14.rel <- merge_phyloseq(ps10.rel, ps14.rel)

obj <- ps10.14.rel
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Generate the matrices you care about.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
head(d50)
#Start cleaning up the data. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#Now will need to break up the strings so that we can identify the samples and make sure we're getting the comparison of interest. 
#library(tidyr)
GU50.m.update <- GU50.m %>%
  separate(Var1, c("cage1", "ear1", "day1"), sep = "_", remove=FALSE)
GU50.m.update$mouse1 <- paste(GU50.m.update$cage1, GU50.m.update$ear1, sep="_")
GU50.m.update <- GU50.m.update %>%
  separate(Var2, c("cage2", "ear2", "day2"), sep = "_", remove=FALSE)
GU50.m.update$mouse2 <- paste(GU50.m.update$cage2, GU50.m.update$ear2, sep="_")

GU50.sd.1 <- GU50.sd[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                       GU50.m.update$day1=="D10" &
                       GU50.m.update$day2=="D14",]
GU50.sd.1 <- GU50.sd.1 %>% drop_na(value)
#
GU50.sd.2 <- GU50.m.update[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                             GU50.m.update$day1=="D14" &
                             GU50.m.update$day2=="D10",]
GU50.sd.2 <- GU50.sd.2 %>% drop_na(value)
#
GU50.sd.2rename <- GU50.sd.2 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1, 
    day1=day2, 
    day2=day1)
GU50.sd.1014 <- rbind(GU50.sd.1, GU50.sd.2rename)
GU50.sd.1014 <- GU50.sd.1014 %>%  #Rename the rownames so that you can merge with metadata df.  
  dplyr::rename(row=Var2)

#
#            ##### from D14 to D28 #####
ps14.28.rel <- merge_phyloseq(ps14.rel, ps28.rel)

obj <- ps14.28.rel
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Generate the matrices you care about.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
head(d50)
#Start cleaning up the data. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#Now will need to break up the strings so that we can identify the samples and make sure we're getting the comparison of interest. 
#library(tidyr)
GU50.m.update <- GU50.m %>%
  separate(Var1, c("cage1", "ear1", "day1"), sep = "_", remove=FALSE)
GU50.m.update$mouse1 <- paste(GU50.m.update$cage1, GU50.m.update$ear1, sep="_")
GU50.m.update <- GU50.m.update %>%
  separate(Var2, c("cage2", "ear2", "day2"), sep = "_", remove=FALSE)
GU50.m.update$mouse2 <- paste(GU50.m.update$cage2, GU50.m.update$ear2, sep="_")

GU50.sd.1 <- GU50.sd[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                       GU50.m.update$day1=="D14" &
                       GU50.m.update$day2=="D28",]
GU50.sd.1 <- GU50.sd.1 %>% drop_na(value)
#
GU50.sd.2 <- GU50.m.update[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                             GU50.m.update$day1=="D28" &
                             GU50.m.update$day2=="D14",]
GU50.sd.2 <- GU50.sd.2 %>% drop_na(value)
#
GU50.sd.2rename <- GU50.sd.2 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1, 
    day1=day2, 
    day2=day1)
GU50.sd.1428 <- rbind(GU50.sd.1, GU50.sd.2rename)
GU50.sd.1428 <- GU50.sd.1428 %>%  #Rename the rownames so that you can merge with metadata df.  
  dplyr::rename(row=Var2)

#
#            ##### from D28 to D41 #####
ps28.41.rel <- merge_phyloseq(ps28.rel, ps41.rel)

obj <- ps28.41.rel
tree <- phyloseq::phy_tree(obj)
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Generate the matrices you care about.
d50 <- gunif[, , "d_0.5"] #Pull the distance matrix. 
head(d50)
#Start cleaning up the data. 
d50[upper.tri(d50)] = NA #This is making sure that there are not duplicate comparisons when we melt the matrix.
GU50.m <- melt(d50, na.rm=TRUE)
GU50.m = GU50.m %>% #Removing all of the 0 values that arise from comparing a sample to itself. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
#Now will need to break up the strings so that we can identify the samples and make sure we're getting the comparison of interest. 
#library(tidyr)
GU50.m.update <- GU50.m %>%
  separate(Var1, c("cage1", "ear1", "day1"), sep = "_", remove=FALSE)
GU50.m.update$mouse1 <- paste(GU50.m.update$cage1, GU50.m.update$ear1, sep="_")
GU50.m.update <- GU50.m.update %>%
  separate(Var2, c("cage2", "ear2", "day2"), sep = "_", remove=FALSE)
GU50.m.update$mouse2 <- paste(GU50.m.update$cage2, GU50.m.update$ear2, sep="_")

GU50.sd.1 <- GU50.sd[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                       GU50.m.update$day1=="D28" &
                       GU50.m.update$day2=="D41",]
GU50.sd.1 <- GU50.sd.1 %>% drop_na(value)
#
GU50.sd.2 <- GU50.m.update[GU50.m.update$mouse1==GU50.m.update$mouse2 & 
                             GU50.m.update$day1=="D41" &
                             GU50.m.update$day2=="D28",]
GU50.sd.2 <- GU50.sd.2 %>% drop_na(value)
#
GU50.sd.2rename <- GU50.sd.2 %>%  #Rename the columns so that all the "pre" samples are 
  dplyr::rename(
    Var1=Var2, 
    Var2=Var1, 
    day1=day2, 
    day2=day1)
GU50.sd.2841 <- rbind(GU50.sd.1, GU50.sd.2rename)
GU50.sd.2841 <- GU50.sd.2841 %>%  #Rename the rownames so that you can merge with metadata df.  
  dplyr::rename(row=Var2)

#

#      ##### Merge beta diversity values together with everything else. #####
#Pull all the GUniFrac values into one df. 
GU50resist <- rbind(GU50.sd.pre1, GU50.sd.110, GU50.sd.1014, GU50.sd.1428, GU50.sd.2841)
GU50resist <- GU50resist %>% rename("GU50" = "value")

#Now merge them with the main resistance df. 
colnames(dfResist)
#Merge the resistance main df and beta div together. 
dfResistwhole <- merge(dfResist, GU50resist, by = "row", all.x = TRUE)
#Replace NAs in beta div pretreatment timepoints with 0 to give the starting point. 
dfResistwhole$GU50[is.na(dfResistwhole$GU50)] <- 0
#

#      ##### Add MDT column for colors #####
dfResist$MDT <- paste(dfResist$microb, dfResist$day_num, dfResist$treatment, sep="_")
dfResist$MDT <- factor(dfResist$MDT, levels=c('LO_-5_worm', 'LO_1_worm', 'LO_10_worm', 'LO_14_worm', 
                                            'LO_28_worm', 'LO_31_worm', 'LO_41_worm', 'HI_-5_worm', 
                                            'HI_1_worm', 'HI_10_worm', 'HI_14_worm', 'HI_28_worm', 
                                            'HI_31_worm', 'HI_41_worm', 'LO_-5_control', 'LO_1_control', 
                                            'LO_10_control', 'LO_14_control', 'LO_28_control', 'LO_31_control', 
                                            'LO_41_control', 'HI_-5_control', 'HI_1_control', 'HI_10_control', 
                                            'HI_14_control', 'HI_28_control', 'HI_31_control', 'HI_41_control'))

#
#   ###### Run stats on resistance using adjacent timepoints as intervals ######
#      ##### Break dataset up into timepoints #####

dfD1resist <- dfResistwhole[dfResistwhole$day_num==1,]
dfD10resist <- dfResistwhole[dfResistwhole$day_num==10,]
dfD14resist <- dfResistwhole[dfResistwhole$day_num==14,]
dfD28resist <- dfResistwhole[dfResistwhole$day_num==28,]
dfD41resist <- dfResistwhole[dfResistwhole$day_num==41,]

#      ##### Example stats testing for resistance #####
#          ##### Richness wilcoxons for all adjacent timepoints.  #####
summary(df1$ChangeASVs[df1$microb == "LO"])
sd(df1$ChangeASVs[df1$microb == "LO"])

#Testing change in ASV counts between pretreatment and D1. 
df1 <- dfD1resist[dfD1resist$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 275, p-value = 0.8003
df2 <- dfD1resist[dfD1resist$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 119.5, p-value = 0.3553
df3 <- dfD1resist[dfD1resist$microb == "HI",]
wilcox.test(ChangeASVs ~ treatment, data=df3)
#W = 22.5, p-value = 0.2682
ps <- c(0.8003 , 0.3553, 0.2682)
p.adjust(ps, "holm")
## 0.8046 0.8046 0.8046

#Testing change in ASV counts between D1 and D10. 
df1 <- dfD10resist[dfD10resist$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 463.5, p-value = 0.0004032
df2 <- dfD10resist[dfD10resist$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 26.5, p-value = 0.005712
df3 <- dfD10resist[dfD10resist$microb == "HI",]
wilcox.test(ChangeASVs ~ treatment, data=df3)
#W = 64, p-value = 0.01917
ps <- c(0.0004032 , 0.005712, 0.01917)
p.adjust(ps, "holm")
## 0.0012096 0.0114240 0.0191700


#Testing change in ASV counts between D10 and D14. 
df1 <- dfD14resist[dfD14resist$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 231, p-value = 0.8461
df2 <- dfD14resist[dfD14resist$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 105, p-value = 0.733
df3 <- dfD14resist[dfD14resist$microb == "HI",]
wilcox.test(ChangeASVs ~ treatment, data=df3)
#W = 16, p-value = 0.1764
ps <- c(0.8461 , 0.733, 0.1764)
p.adjust(ps, "holm")
## 1.0000 1.0000 0.5292

#Testing change in ASV counts between D14 and D28. 
df1 <- dfD28resist[dfD28resist$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 58, p-value = 0.005894
df2 <- dfD28resist[dfD28resist$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 30.5, p-value = 0.2748
#
ps <- c(0.005894 , 0.2748)
p.adjust(ps, "holm")
## 0.011788, 0.274800

#Testing change in evenness between D28 and D41
df1 <- dfD41resist[dfD41resist$treatment == "worm",]
wilcox.test(ChangeASVs ~ microb, data=df1)
#W = 7, p-value = 0.008803
df2 <- dfD41resist[dfD41resist$microb == "LO",]
wilcox.test(ChangeASVs ~ treatment, data=df2)
#W = 73, p-value = 0.032
ps <- c(0.008803 , 0.032)
p.adjust(ps, "holm")
## 0.017606, 0.032000


#          ##### GUniFrac wilcoxons for all adjacent timepoints.  #####

#Testing change in GUniFrac distance between pretreatment and D1. 
df1 <- dfD1resist[dfD1resist$treatment == "worm",]
wilcox.test(GU50 ~ microb, data=df1)
#W = 148, p-value = 0.00408
df2 <- dfD1resist[dfD1resist$microb == "LO",]
wilcox.test(GU50 ~ treatment, data=df2)
#W = 96, p-value = 1
df3 <- dfD1resist[dfD1resist$microb == "HI",]
wilcox.test(GU50 ~ treatment, data=df3)
#W = 41, p-value = 0.712
ps <- c(0.00408 , 1, 0.712)
p.adjust(ps, "holm")
## 0.01224 1.00000 1.00000

#Testing change in GUniFrac distance between D1 and D10. 
df1 <- dfD10resist[dfD10resist$treatment == "worm",]
wilcox.test(GU50 ~ microb, data=df1)
#W = 320, p-value = 0.528
df2 <- dfD10resist[dfD10resist$microb == "LO",]
wilcox.test(GU50 ~ treatment, data=df2)
#W = 59, p-value = 0.1473
df3 <- dfD10resist[dfD10resist$microb == "HI",]
wilcox.test(GU50 ~ treatment, data=df3)
#W = 56, p-value = 0.09815
ps <- c(0.528 , 0.1473, 0.09815)
p.adjust(ps, "holm")
## 0.52800 0.29460 0.29445


#Testing change in GUniFrac distance between D10 and D14. 
df1 <- dfD14resist[dfD14resist$treatment == "worm",]
wilcox.test(GU50 ~ microb, data=df1)
#W = 200, p-value = 0.3715
df2 <- dfD14resist[dfD14resist$microb == "LO",]
wilcox.test(GU50 ~ treatment, data=df2)
#W = 58, p-value = 0.1361
df3 <- dfD14resist[dfD14resist$microb == "HI",]
wilcox.test(GU50 ~ treatment, data=df3)
#W = 34, p-value = 0.7363
ps <- c(0.3715 , 0.1361, 0.7363)
p.adjust(ps, "holm")
## 0.7430 0.4083 0.7430

#Testing change in GUniFrac distance between D14 and D28. 
df1 <- dfD28resist[dfD28resist$treatment == "worm",]
wilcox.test(GU50 ~ microb, data=df1)
#W = 34, p-value = 0.7363
df2 <- dfD28resist[dfD28resist$microb == "LO",]
wilcox.test(GU50 ~ treatment, data=df2)
#W = 22, p-value = 0.07965
#
ps <- c(0.7363 , 0.07965)
p.adjust(ps, "holm")
## 0.7363 0.1593

#Testing change in evenness between D28 and D41
df1 <- dfD41resist[dfD41resist$treatment == "worm",]
wilcox.test(GU50 ~ microb, data=df1)
#W = 23, p-value = 0.23
df2 <- dfD41resist[dfD41resist$microb == "LO",]
wilcox.test(GU50 ~ treatment, data=df2)
#W = 33, p-value = 0.3809
ps <- c(0.23 , 0.3809)
p.adjust(ps, "holm")
## 0.46 0.46



#   ##### Example boxplots of resistance intervals #####
#Note that decided it doesn't make sense to make line plots of adjacent timepoints. Just stick with boxplots. 

#      ##### ASV Richness #####
#Boxplot change in richness pretreatment to D1
pRichnesspre1 <- ggplot(dfD1resist, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  #ylim(min(dfD1resist$ChangeASVs), max(dfD1resist$ChangeASVs)+0.1*max(dfD1resist$ChangeASVs)) +
  ylim(-100, 120) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pretreatment-D1") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichnesspre1
#
#Plotting change in richness D1 to D10
pRichness110 <- ggplot(dfD10resist, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  #ylim(min(dfD10resist$ChangeASVs), max(dfD10resist$ChangeASVs)+0.1*max(dfD10resist$ChangeASVs)) +
  ylim(-100, 120) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness D1-D10") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichness110
#
#Plotting change in richness D10 to D14
pRichness1014 <- ggplot(dfD14resist, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  #ylim(min(dfD14resist$ChangeASVs), max(dfD14resist$ChangeASVs)+0.1*max(dfD14resist$ChangeASVs)) +
  ylim(-100, 120) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness D10-D14") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichness1014
#
#Boxplot change in richness D14 to D28
pRichness1428  <- ggplot(dfD28resist, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964")) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD28resist$ChangeASVs), max(dfD28resist$ChangeASVs)+0.1*max(dfD28resist$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness D14-D28") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichness1428
#
#Boxplot change in richness D28 to D41
pRichness2841  <- ggplot(dfD41resist, aes(x=MT, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964")) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfResil2841$ChangeASVs), max(dfResil2841$ChangeASVs)+0.1*max(dfResil2841$ChangeASVs)) +
#  ylim(min(dfD41resist$ChangeASVs), max(dfD41resist$ChangeASVs)+0.1*max(dfD41resist$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness D28-D41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pRichness2841
#

#         ##### Save the plots for Richness #####
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_ChangeRichnesspre1_boxplot_woD31.pdf", pRichnesspre1, height = 5, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_ChangeRichness110_boxplot_woD31.pdf", pRichness110, height = 5, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_ChangeRichness1014_boxplot_woD31.pdf", pRichness1014, height = 5, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_ChangeRichness1428_boxplot_woD31.pdf", pRichness1428, height = 5, width = 3, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_ChangeRichness2841_boxplot_woD31.pdf", pRichness2841, height = 5, width = 3, units = "in")
#
#
#      ##### GUniFrac #####
#Boxplot change in richness pretreatment to D1
pGUniFracpre1 <- ggplot(dfD1resist, aes(x=MT, y=GU50)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  #ylim(min(dfD1resist$GU50), max(dfD1resist$GU50)+0.1*max(dfD1resist$GU50)) +
  ylim(0.05, 0.7)+
  ylab("Change in generalized UniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac-Curtis distance pretreatment-D1") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFracpre1
#
#Plotting change in richness D1 to D10
pGUniFrac110 <- ggplot(dfD10resist, aes(x=MT, y=GU50)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  #ylim(min(dfD10resist$GU50), max(dfD10resist$GU50)+0.1*max(dfD10resist$GU50)) +
  ylim(0.05, 0.7)+
  ylab("Change in generalized UniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac-Curtis distance D1-D10") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFrac110
#
#Plotting change in richness D10 to D14
pGUniFrac1014 <- ggplot(dfD14resist, aes(x=MT, y=GU50)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964","#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964","#FFFFFF")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  #ylim(min(dfD14resist$GU50), max(dfD14resist$GU50)+0.1*max(dfD14resist$GU50)) +
  ylim(0.05, 0.7)+
  ylab("Change in generalized UniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac-Curtis distance D10-D14") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFrac1014
#
#Boxplot change in richness D14 to D28
pGUniFrac1428  <- ggplot(dfD28resist, aes(x=MT, y=GU50)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964")) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD28resist$GU50), max(dfD28resist$GU50)+0.1*max(dfD28resist$GU50)) +
  #ylim(0.05, 0.5)+
  ylab("Change in generalized UniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac-Curtis distance D14-D28") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFrac1428
#
#Boxplot change in richness D28 to D41
pGUniFrac2841  <- ggplot(dfD41resist, aes(x=MT, y=GU50)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964")) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD41resist$GU50), max(dfResil2841$ChangeGUniFrac)+0.1*max(dfResil2841$ChangeGUniFrac)) +
#  ylim(min(dfD41resist$GU50), max(dfD41resist$GU50)+0.1*max(dfD41resist$GU50)) +
  #ylim(0.05, 0.5)+
  ylab("Change in generalized UniFrac distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac-Curtis distance D28-D41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) 
pGUniFrac2841
#

#         ##### Save the plots for GUniFrac-Curtis distance #####
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resistance/psexp3000_GUniFracpre1_boxplot_woD31.pdf", pGUniFracpre1, height = 5, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resistance/psexp3000_GUniFrac110_boxplot_woD31.pdf", pGUniFrac110, height = 5, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resistance/psexp3000_GUniFrac1014_boxplot_woD31.pdf", pGUniFrac1014, height = 5, width = 4, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resistance/psexp3000_GUniFrac1428_boxplot_woD31.pdf", pGUniFrac1428, height = 5, width = 3, units = "in")
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resistance/psexp3000_GUniFrac2841_boxplot_woD31.pdf", pGUniFrac2841, height = 5, width = 3, units = "in")
#
#
######## Return to RAW numbers from all timepoints into the same type of df -- RAW NUMBERS #####
dfexp3000no31

library(reshape2)
#Melt for ASV counts. 
dfmeltobs <- melt(dfmerge, id.vars=c(1,6:10), measure.vars=c(17, 24, 28, 32, 36, 40), 
                  value.name="RawASVs", variable.name="time")
View(dfmeltobs)

#Melt evenness counts
dfmeltevenness <- melt(dfmerge, id.vars=c(1,6:10), measure.vars=c(23, 27, 31, 35, 39, 43), 
                       value.name="RawEvenness", variable.name="time")

#Stick the change in shannon and faiths on the end of the obs to create back to 1 df. 
dfraw <- cbind(dfmeltobs, dfmeltevenness$RawEvenness)
dfraw <- dfraw %>% dplyr::rename(RawEvenness='dfmeltevenness$RawEvenness')
dfraw$MT <- paste(dfraw$microb, dfraw$treatment, sep="_")
dfraw$MT <- factor(dfraw$MT, levels=c("LO_worm", "LO_control", "HI_worm", "HI_control"))
dfraw$day_num <- as.character(dfraw$time)
dfraw <- dfraw %>% tidyr::drop_na(RawASVs)

#Create a numeric day column. 
dfraw <- dfraw %>% 
  mutate(day_num = replace(day_num, day_num == 'obs.pre', -5)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs.D1', 1)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs.D10', 10)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs.D14', 14)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs.D28', 28)) %>%
  mutate(day_num = replace(day_num, day_num == 'obs.D41', 41))
#Now make the column numeric. 
dfraw$day_num <- as.numeric(dfraw$day_num)

View(dfraw) #Make sure it all came out ok. 
#
#Get a quick view of the dataset
ggplot(dfraw, aes(x=day_num, y=RawASVs)) + geom_point(aes(color=MT)) + theme_bw()

#   ##### Write out csv of the df for safekeeping.  #####

write.csv(dfraw, "AlphaDiv/RawAlphaDiv_noD31.csv")
#
#   ##### Remove incomplete data points and really huge outliers from data.  #####
#Remove the rows for D28 for all HI and remove HI sham D31 bc only 1 point. Lacking data.  
dfraw1 <- dfraw[!(dfraw$day_num == "28" & dfraw$MT == "HI_control"),]
dfraw2 <- dfraw1[!(dfraw1$day_num == "41" & dfraw1$mouse == "21_LP"),]

#Run anovas on the days of interest. 
dfraw2pre <- dfraw2[dfraw2$day_num==-5,]
dfraw2D1 <- dfraw2[dfraw2$day_num==1,]
dfraw2D10 <- dfraw2[dfraw2$day_num==10,]
dfraw2D14 <- dfraw2[dfraw2$day_num==14,]
dfraw2D28 <- dfraw2[dfraw2$day_num==28,]
dfraw2D31 <- dfraw2[dfraw2$day_num==31,]
dfraw2D41 <- dfraw2[dfraw2$day_num==41,]
#
#Testing for normality
hist(dfraw2D10$c)
shapiro.test(dfraw2D10$RawASVs)
#Values are highly significant. Should not use parametrics

#            ##### Run wilcoxons on each timepoint of interest for raw richness. ######
#Just picked the couple of timepoints of interest during infection. 
#Testing change pretreatment to D10. 
df1 <- dfraw2D10[dfraw2D10$treatment == "worm",]
wilcox.test(RawASVs ~ microb, data=df1)
#W = 0, p-value = 6.178e-09 
df2 <- dfraw2D10[dfraw2D10$microb == "LO",]
wilcox.test(RawASVs ~ treatment, data=df2)
#W = 53, p-value = 0.08872
df3 <- dfraw2D10[dfraw2D10$microb == "HI",]
wilcox.test(RawASVs ~ treatment, data=df3)
#W = 68, p-value = 0.007314
#
ps <- c(0.000000006178, 0.08872, 0.007314)
p.adjust(ps, "holm")
##1.8534e-08, 8.8720e-02, 1.4628e-02


#Testing change pretreatment to D14. 
df1 <- dfraw2D14[dfraw2D14$treatment == "worm",]
wilcox.test(RawASVs ~ microb, data=df1)
#W = 0, p-value = 4.474e-08
df2 <- dfraw2D14[dfraw2D14$microb == "LO",]
wilcox.test(RawASVs ~ treatment, data=df2)
#W = 55.5, p-value = 0.1088
df3 <- dfraw2D14[dfraw2D14$microb == "HI",]
wilcox.test(RawASVs ~ treatment, data=df3)
#W = 42, p-value = 0.2483
#
ps <- c(0.00000004474, 0.1088, 0.2483)
p.adjust(ps, "holm")
## 1.3422e-07, 2.1760e-01, 2.4830e-01

#      ##### Create RAW means and sd for plotting and comparison of ASVs #####

#Take a quick look at the data to make sure stuff looks like it's supposed to. 
ggplot(dfraw2, aes(x=day_num, y=RawASVs)) +
  geom_point(aes(colour=MT)) +
  theme_bw()


#Ok, now start the subsetting and calculations to make your moving averages plot. 
dfhitreat <- dfraw2[dfraw2$microb=="HI" & dfraw2$treatment=="worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(RawASVs)
dfhictrl <- dfraw2[dfraw2$microb=="HI" & dfraw2$treatment=="control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(RawASVs)
dflotreat <- dfraw2[dfraw2$microb=="LO" & dfraw2$treatment=="worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(RawASVs)
dfloctrl <- dfraw2[dfraw2$microb=="LO" & dfraw2$treatment=="control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(RawASVs)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawASVs = mean(RawASVs, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawASVs = mean(RawASVs, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawASVs = mean(RawASVs, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawASVs = mean(RawASVs, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawASVs, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawASVs, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawASVs, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawASVs, na.rm=TRUE))

#
dfhitreat.mean$SD <- dfhitreat.sd$SD
dfhictrl.mean$SD <- dfhictrl.sd$SD
dflotreat.mean$SD <- dflotreat.sd$SD
dfloctrl.mean$SD <- dfloctrl.sd$SD


# add sample size and calculate SEM 
hitreat.N <- as.data.frame(plyr::count(dfhitreat$day_num)) 
hictrl.N <- as.data.frame(plyr::count(dfhictrl$day_num)) 
lotreat.N <- as.data.frame(plyr::count(dflotreat$day_num)) 
loctrl.N <- as.data.frame(plyr::count(dfloctrl$day_num)) 

# add N to dataframes 
dfhitreat.mean$N <- hitreat.N$freq
dfhictrl.mean$N <- hictrl.N$freq
dflotreat.mean$N <- lotreat.N$freq
dfloctrl.mean$N <- loctrl.N$freq

# add SE of each group to the mean. 

dfhitreat.mean$SEH <- dfhitreat.mean$meanRawASVs + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanRawASVs + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanRawASVs + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanRawASVs + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanRawASVs - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanRawASVs - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanRawASVs - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanRawASVs - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as asv-specific dfs. 
dfhitreat.mean.asv <- dfhitreat.mean
dfhictrl.mean.asv <- dfhictrl.mean 
dflotreat.mean.asv <- dflotreat.mean
dfloctrl.mean.asv <- dfloctrl.mean 

#         ##### Plot RAW RICHNESS as average over time with error bars #####
ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.asv$meanRawASVs, x=dfhitreat.mean.asv$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.asv$meanRawASVs, x=dfhictrl.mean.asv$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.asv$SEL, ymax=dfhitreat.mean.asv$SEH, y=dfhitreat.mean.asv$meanRawASVs, x=dfhitreat.mean.asv$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.asv$SEL, ymax=dfhictrl.mean.asv$SEH, y=dfhictrl.mean.asv$meanRawASVs, x=dfhictrl.mean.asv$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.asv$meanRawASVs, x=dfhitreat.mean.asv$day_num), shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.asv$meanRawASVs, x=dfhictrl.mean.asv$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.asv$meanRawASVs, x=dflotreat.mean.asv$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.asv$meanRawASVs, x=dfloctrl.mean.asv$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.asv$SEL, ymax=dflotreat.mean.asv$SEH, y=dflotreat.mean.asv$meanRawASVs, x=dflotreat.mean.asv$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.asv$SEL, ymax=dfloctrl.mean.asv$SEH, y=dfloctrl.mean.asv$meanRawASVs, x=dfloctrl.mean.asv$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.asv$meanRawASVs, x=dflotreat.mean.asv$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.asv$meanRawASVs, x=dfloctrl.mean.asv$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("ASV richness") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average ASV richness across all timepoints") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

ggsave("AlphaDiv/psexp3000_ASVRichness_LineplotwAverages_woD31.pdf", units="in", height=5, width=6)
#



#         ##### Plot as average over time with error bars and colors at different timepoints #####
locols = c('#d6a3f0' , '#c275eb' , '#A11FE5', '#7E05BE' , '#5C008D', '#300049')

hicols = c("#f5e1a6" , "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" ,"#644C00")

library(ggnewscale) #To be able to set multiple color scales in a single plot

RawASVChangeplot <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + 
  geom_line(aes(y=dfhitreat.mean.asv$meanRawASVs, x=dfhitreat.mean.asv$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.asv$meanRawASVs, x=dfhictrl.mean.asv$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.asv$SEL, ymax=dfhitreat.mean.asv$SEH, y=dfhitreat.mean.asv$meanRawASVs, x=dfhitreat.mean.asv$day_num, color=as.character(dfhitreat.mean.asv$day_num)), alpha=1,  width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.asv$SEL, ymax=dfhictrl.mean.asv$SEH, y=dfhictrl.mean.asv$meanRawASVs, x=dfhictrl.mean.asv$day_num, color=as.character(dfhictrl.mean.asv$day_num)), alpha=1, width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.asv$meanRawASVs, x=dfhitreat.mean.asv$day_num, color=as.character(dfhitreat.mean.asv$day_num)),shape=15, alpha=1, size=3, na.rm = TRUE) + 
  geom_point(aes(y=dfhictrl.mean.asv$meanRawASVs, x=dfhictrl.mean.asv$day_num, color=as.character(dfhictrl.mean.asv$day_num)), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE) + 
  scale_color_manual(values=hicols) +
  new_scale_color() +
  geom_line(aes(y=dflotreat.mean.asv$meanRawASVs, x=dflotreat.mean.asv$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.asv$meanRawASVs, x=dfloctrl.mean.asv$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.asv$SEL, ymax=dflotreat.mean.asv$SEH, y=dflotreat.mean.asv$meanRawASVs, x=dflotreat.mean.asv$day_num, color=as.character(dflotreat.mean.asv$day_num)), alpha=1,  width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.asv$SEL, ymax=dfloctrl.mean.asv$SEH, y=dfloctrl.mean.asv$meanRawASVs, x=dfloctrl.mean.asv$day_num, color=as.character(dfloctrl.mean.asv$day_num)), alpha=1, width=0.8) + 
  geom_point(aes(y=dflotreat.mean.asv$meanRawASVs, x=dflotreat.mean.asv$day_num, color=as.character(dflotreat.mean.asv$day_num)),shape=16, alpha=1, size=3, na.rm = TRUE) + 
  geom_point(aes(y=dfloctrl.mean.asv$meanRawASVs, x=dfloctrl.mean.asv$day_num, color=as.character(dfloctrl.mean.asv$day_num)), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE) + 
  scale_color_manual(values=locols) +
  ylab("ASV richness") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average ASV") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

RawASVChangeplot

ggsave("AlphaDiv/psexp3000_ASVRichness_LineplotwAverages_allcolors_woD31.pdf", RawASVChangeplot, units="in", height=5, width=6)
#
#      ##### Create means and sd for plotting and comparison of Evenness #####

#Take a quick look at the data to make sure stuff looks like it's supposed to. 
ggplot(dfraw2, aes(x=day_num, y=RawEvenness)) +
  geom_point(aes(colour=MT)) +
  theme_bw()


#Ok, now start the subsetting and calculations to make your moving averages plot. 
dfhitreat <- dfraw2[dfraw2$microb=="HI" & dfraw2$treatment=="worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(RawEvenness)
dfhictrl <- dfraw2[dfraw2$microb=="HI" & dfraw2$treatment=="control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(RawEvenness)
dflotreat <- dfraw2[dfraw2$microb=="LO" & dfraw2$treatment=="worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(RawEvenness)
dfloctrl <- dfraw2[dfraw2$microb=="LO" & dfraw2$treatment=="control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(RawEvenness)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawEvenness = mean(RawEvenness, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawEvenness = mean(RawEvenness, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawEvenness = mean(RawEvenness, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanRawEvenness = mean(RawEvenness, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawEvenness, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawEvenness, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawEvenness, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(RawEvenness, na.rm=TRUE))

#
dfhitreat.mean$SD <- dfhitreat.sd$SD
dfhictrl.mean$SD <- dfhictrl.sd$SD
dflotreat.mean$SD <- dflotreat.sd$SD
dfloctrl.mean$SD <- dfloctrl.sd$SD


# add sample size and calculate SEM 
hitreat.N <- as.data.frame(plyr::count(dfhitreat$day_num)) 
hictrl.N <- as.data.frame(plyr::count(dfhictrl$day_num)) 
lotreat.N <- as.data.frame(plyr::count(dflotreat$day_num)) 
loctrl.N <- as.data.frame(plyr::count(dfloctrl$day_num)) 

# add N to dataframes 
dfhitreat.mean$N <- hitreat.N$freq
dfhictrl.mean$N <- hictrl.N$freq
dflotreat.mean$N <- lotreat.N$freq
dfloctrl.mean$N <- loctrl.N$freq

# add SE of each group to the mean. 

dfhitreat.mean$SEH <- dfhitreat.mean$meanRawEvenness + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanRawEvenness + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanRawEvenness + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanRawEvenness + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanRawEvenness - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanRawEvenness - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanRawEvenness - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanRawEvenness - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as asv-specific dfs. 
dfhitreat.mean.asv <- dfhitreat.mean
dfhictrl.mean.asv <- dfhictrl.mean 
dflotreat.mean.asv <- dflotreat.mean
dfloctrl.mean.asv <- dfloctrl.mean 

#         ##### Plot as average over time with error bars #####
ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.asv$meanRawEvenness, x=dfhitreat.mean.asv$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.asv$meanRawEvenness, x=dfhictrl.mean.asv$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.asv$SEL, ymax=dfhitreat.mean.asv$SEH, y=dfhitreat.mean.asv$meanRawEvenness, x=dfhitreat.mean.asv$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.asv$SEL, ymax=dfhictrl.mean.asv$SEH, y=dfhictrl.mean.asv$meanRawEvenness, x=dfhictrl.mean.asv$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.asv$meanRawEvenness, x=dfhitreat.mean.asv$day_num), shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.asv$meanRawEvenness, x=dfhictrl.mean.asv$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.asv$meanRawEvenness, x=dflotreat.mean.asv$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.asv$meanRawEvenness, x=dfloctrl.mean.asv$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.asv$SEL, ymax=dflotreat.mean.asv$SEH, y=dflotreat.mean.asv$meanRawEvenness, x=dflotreat.mean.asv$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.asv$SEL, ymax=dfloctrl.mean.asv$SEH, y=dfloctrl.mean.asv$meanRawEvenness, x=dfloctrl.mean.asv$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.asv$meanRawEvenness, x=dflotreat.mean.asv$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.asv$meanRawEvenness, x=dfloctrl.mean.asv$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("Pielou's evenness") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average Evenness across all timepoints") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

ggsave("AlphaDiv/psexp3000_Evenness_LineplotwAverages_woD31.pdf", units="in", height=5, width=6)
#






###############