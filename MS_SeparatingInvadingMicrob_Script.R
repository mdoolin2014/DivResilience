#Looking for invaders into microbiotas during each sampling interval.  

psexp3000no31 <- readRDS("~/Documents/1U4U_16S/GG2_reclassification/R_objects/psexp3000no31.rds")
#

##### FINDING ASVs THAT INVADED AT ADJACENT TIMEPOINTS #####
#Have already created the phyloseq objects needed (individual timepoints by each subgroup) in Start script. 

#   ##### For pretreatment to D1 #####

lotreatearly <- psprelotreat
loctrlearly <- pspreloctrl
hitreatearly <- psprehitreat
hictrlearly <- psprehictrl

lotreattax <- as.list(rownames(otu_table(lotreatearly)))
loctrltax <- as.list(rownames(otu_table(loctrlearly)))
hitreattax <- as.list(rownames(otu_table(hitreatearly)))
hictrltax <- as.list(rownames(otu_table(hictrlearly)))

lotreatlate <- ps1lotreat
loctrllate <- ps1loctrl
hitreatlate <- ps1hitreat
hictrllate <- ps1hictrl

length(loctrltax) #105 ASVs
tmpa <- subset(otu_table(lotreatlate), !rownames(otu_table(lotreatlate)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatlate), metadata=sample_data(lotreatlate))
#25 taxa remaining

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrllate), !rownames(otu_table(loctrllate)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrllate), metadata=sample_data(loctrllate))
#13 taxa remaining

tmpc <- subset(otu_table(hitreatlate), !rownames(otu_table(hitreatlate)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatlate), metadata=sample_data(hitreatlate))
#45 taxa remaining

tmpd <- subset(otu_table(hictrllate), !rownames(otu_table(hictrllate)) %in% hictrltax)
tmpd1 <- merge_phyloseq(tmpd, taxonomy=tax_table(hictrllate), metadata=sample_data(hictrllate))
#53 taxa remaining

#Now let's bring everything together into one large object. 
pspre1invaders <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tmpd1, tree = phy_tree(psexp3000no31))
phyobj <- pspre1invaders
sample_data(phyobj)$comparison <- "pre1"
tmp <- data.frame(estimate_richness(phyobj, measures="Observed"))
sample_data(phyobj)$richnessInvaders <- tmp$Observed
sample_data(phyobj)$readsInvaders <- sample_sums(phyobj)
sample_data(phyobj)$sample <- rownames(sample_data(phyobj))
pspre1invaders <- phyobj
#View(sample_data(pspre1invaders))
#

#Now pull out the tax_tables and label them based on which groups they were invading. 
#Getting this ready to merge with the Differential Abundance tables for updated plotting.  
tmptaxa <- cbind(data.frame(tax_table(tmpa1), rowSums(otu_table(tmpa1))), "LO_worm_pre1")
colnames(tmptaxa)[8] <- "ReadSums"
colnames(tmptaxa)[9] <- "MTD"
tmptaxb <- cbind(data.frame(tax_table(tmpb1), rowSums(otu_table(tmpb1))), "LO_sham_pre1")
colnames(tmptaxb)[8] <- "ReadSums"
colnames(tmptaxb)[9] <- "MTD"
tmptaxc <- cbind(data.frame(tax_table(tmpc1), rowSums(otu_table(tmpc1))), "HI_worm_pre1")
colnames(tmptaxc)[8] <- "ReadSums"
colnames(tmptaxc)[9] <- "MTD"
tmptaxd <- cbind(data.frame(tax_table(tmpd1), rowSums(otu_table(tmpd1))), "HI_sham_pre1")
colnames(tmptaxd)[8] <- "ReadSums"
colnames(tmptaxd)[9] <- "MTD"

taxpre1 <- rbind(tmptaxa, tmptaxb, tmptaxc, tmptaxd)

#
#
#   ##### For D1 to D10 #####

lotreatearly <- ps1lotreat
loctrlearly <- ps1loctrl
hitreatearly <- ps1hitreat
hictrlearly <- ps1hictrl

lotreattax <- as.list(rownames(otu_table(lotreatearly)))
loctrltax <- as.list(rownames(otu_table(loctrlearly)))
hitreattax <- as.list(rownames(otu_table(hitreatearly)))
hictrltax <- as.list(rownames(otu_table(hictrlearly)))

lotreatlate <- ps10lotreat
loctrllate <- ps10loctrl
hitreatlate <- ps10hitreat
hictrllate <- ps10hictrl

length(loctrltax) #105 ASVs
tmpa <- subset(otu_table(lotreatlate), !rownames(otu_table(lotreatlate)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatlate), metadata=sample_data(lotreatlate))
#126 taxa remaining
richa <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "LO_worm_110",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpa1)), 
                         "Family" = data.frame(tax_table(tmpa1))$Family)
richa1 <- merge(richa, tmpdat, by="ASV", all=FALSE)

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrllate), !rownames(otu_table(loctrllate)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrllate), metadata=sample_data(loctrllate))
#18 taxa remaining

tmpc <- subset(otu_table(hitreatlate), !rownames(otu_table(hitreatlate)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatlate), metadata=sample_data(hitreatlate))
#56 taxa remaining
richc <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "HI_worm_110",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpc1)), 
                     "Family" = data.frame(tax_table(tmpc1))$Family)
richc1 <- merge(richc, tmpdat, by="ASV", all=FALSE)
#0 ASVs that invaded were differentially abundant. 

tmpd <- subset(otu_table(hictrllate), !rownames(otu_table(hictrllate)) %in% hictrltax)
tmpd1 <- merge_phyloseq(tmpd, taxonomy=tax_table(hictrllate), metadata=sample_data(hictrllate))
#109 taxa remaining
richd <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "HI_sham_110",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpd1)), 
                     "Family" = data.frame(tax_table(tmpd1))$Family)
richd1 <- merge(richd, tmpdat, by="ASV", all=FALSE)
#1 ASV that invaded were differentially abundant. -- Oscillospiraceae 7f9cf...



#Now let's bring everything together into one large object. 
ps110invaders <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tmpd1, tree = phy_tree(psexp3000no31))
phyobj <- ps110invaders
sample_data(phyobj)$comparison <- "110"
sample_data(phyobj)$readsInvaders <- sample_sums(phyobj)
tmp <- data.frame(estimate_richness(phyobj, measures="Observed"))
sample_data(phyobj)$richnessInvaders <- tmp$Observed
sample_data(phyobj)$sample <- rownames(sample_data(phyobj))
ps110invaders <- phyobj

#View(sample_data(ps110invaders))

#Now extract the taxa from the invaders and quantify, so that I can merge with the diff abund data. 
tmptaxa <- cbind(data.frame(tax_table(tmpa1), rowSums(otu_table(tmpa1))), "LO_worm_110")
colnames(tmptaxa)[8] <- "ReadSums"
colnames(tmptaxa)[9] <- "MTD"
tmptaxb <- cbind(data.frame(tax_table(tmpb1), rowSums(otu_table(tmpb1))), "LO_sham_110")
colnames(tmptaxb)[8] <- "ReadSums"
colnames(tmptaxb)[9] <- "MTD"
tmptaxc <- cbind(data.frame(tax_table(tmpc1), rowSums(otu_table(tmpc1))), "HI_worm_110")
colnames(tmptaxc)[8] <- "ReadSums"
colnames(tmptaxc)[9] <- "MTD"
tmptaxd <- cbind(data.frame(tax_table(tmpd1), rowSums(otu_table(tmpd1))), "HI_sham_110")
colnames(tmptaxd)[8] <- "ReadSums"
colnames(tmptaxd)[9] <- "MTD"

tax110 <- rbind(tmptaxa, tmptaxb, tmptaxc, tmptaxd)
#
#
#   ##### For D10 to D14 #####
lotreatearly <- ps10lotreat
loctrlearly <- ps10loctrl
hitreatearly <- ps10hitreat
hictrlearly <- ps10hictrl

lotreattax <- as.list(rownames(otu_table(lotreatearly)))
loctrltax <- as.list(rownames(otu_table(loctrlearly)))
hitreattax <- as.list(rownames(otu_table(hitreatearly)))
hictrltax <- as.list(rownames(otu_table(hictrlearly)))

lotreatlate <- ps14lotreat
loctrllate <- ps14loctrl
hitreatlate <- ps14hitreat
hictrllate <- ps14hictrl

length(loctrltax) #105 ASVs
tmpa <- subset(otu_table(lotreatlate), !rownames(otu_table(lotreatlate)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatlate), metadata=sample_data(lotreatlate))
#14 taxa remaining

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrllate), !rownames(otu_table(loctrllate)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrllate), metadata=sample_data(loctrllate))
#12 taxa remaining
richb <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "LO_sham_1014",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpb1)), 
                     "Family" = data.frame(tax_table(tmpb1))$Family)
richb1 <- merge(richb, tmpdat, by="ASV", all=FALSE)


tmpc <- subset(otu_table(hitreatlate), !rownames(otu_table(hitreatlate)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatlate), metadata=sample_data(hitreatlate))
#67 taxa remaining
richc <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "HI_worm_1014",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpc1)), 
                     "Family" = data.frame(tax_table(tmpc1))$Family)
richc1 <- merge(richc, tmpdat, by="ASV", all=FALSE)


tmpd <- subset(otu_table(hictrllate), !rownames(otu_table(hictrllate)) %in% hictrltax)
tmpd1 <- merge_phyloseq(tmpd, taxonomy=tax_table(hictrllate), metadata=sample_data(hictrllate))
#54 taxa remaining
richd <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "HI_sham_1014",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpd1)), 
                     "Family" = data.frame(tax_table(tmpd1))$Family)
richd1 <- merge(richd, tmpdat, by="ASV", all=FALSE)


#Now let's bring everything together into one large object. 
ps1014invaders <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tmpd1, tree = phy_tree(psexp3000))
phyobj <- ps1014invaders
sample_data(phyobj)$comparison <- "1014"
sample_data(phyobj)$readsInvaders <- sample_sums(phyobj)
tmp <- data.frame(estimate_richness(phyobj, measures="Observed"))
sample_data(phyobj)$richnessInvaders <- tmp$Observed
sample_data(phyobj)$sample <- rownames(sample_data(phyobj))
ps1014invaders <- phyobj
#View(sample_data(ps1014invaders))


#Now extract the taxa from the invaders and quantify, so that I can merge with the 
tmptaxa <- cbind(data.frame(tax_table(tmpa1), rowSums(otu_table(tmpa1))), "LO_worm_1014")
colnames(tmptaxa)[8] <- "ReadSums"
colnames(tmptaxa)[9] <- "MTD"
tmptaxb <- cbind(data.frame(tax_table(tmpb1), rowSums(otu_table(tmpb1))), "LO_sham_1014")
colnames(tmptaxb)[8] <- "ReadSums"
colnames(tmptaxb)[9] <- "MTD"
tmptaxc <- cbind(data.frame(tax_table(tmpc1), rowSums(otu_table(tmpc1))), "HI_worm_1014")
colnames(tmptaxc)[8] <- "ReadSums"
colnames(tmptaxc)[9] <- "MTD"
tmptaxd <- cbind(data.frame(tax_table(tmpd1), rowSums(otu_table(tmpd1))), "HI_sham_1014")
colnames(tmptaxd)[8] <- "ReadSums"
colnames(tmptaxd)[9] <- "MTD"

tax1014 <- rbind(tmptaxa, tmptaxb, tmptaxc, tmptaxd)
#


#   ##### For D14 to D28 #####

lotreatearly <- ps14lotreat
loctrlearly <- ps14loctrl
hitreatearly <- ps14hitreat

lotreattax <- as.list(rownames(otu_table(lotreatearly)))
loctrltax <- as.list(rownames(otu_table(loctrlearly)))
hitreattax <- as.list(rownames(otu_table(hitreatearly)))

lotreatlate <- ps28lotreat
loctrllate <- ps28loctrl
hitreatlate <- ps28hitreat

length(lotreattax) #
tmpa <- subset(otu_table(lotreatlate), !rownames(otu_table(lotreatlate)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatlate), metadata=sample_data(lotreatlate))
#9 taxa remaining

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrllate), !rownames(otu_table(loctrllate)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrllate), metadata=sample_data(loctrllate))
#17 taxa remaining

tmpc <- subset(otu_table(hitreatlate), !rownames(otu_table(hitreatlate)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatlate), metadata=sample_data(hitreatlate))
#7 taxa remaining

#Now let's bring everything together into one large object. 
ps1428invaders <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tree = phy_tree(psexp3000))
phyobj <- ps1428invaders
sample_data(phyobj)$comparison <- "1428"
sample_data(phyobj)$readsInvaders <- sample_sums(phyobj)
tmp <- data.frame(estimate_richness(phyobj, measures="Observed"))
sample_data(phyobj)$richnessInvaders <- tmp$Observed
sample_data(phyobj)$sample <- rownames(sample_data(phyobj))
ps1428invaders <- phyobj
#View(sample_data(ps1428invaders))

#Now extract the taxa from the invaders and quantify, so that I can merge with the 
tmptaxa <- cbind(data.frame(tax_table(tmpa1), rowSums(otu_table(tmpa1))), "LO_worm_1428")
colnames(tmptaxa)[8] <- "ReadSums"
colnames(tmptaxa)[9] <- "MTD"
tmptaxb <- cbind(data.frame(tax_table(tmpb1), rowSums(otu_table(tmpb1))), "LO_sham_1428")
colnames(tmptaxb)[8] <- "ReadSums"
colnames(tmptaxb)[9] <- "MTD"
tmptaxc <- cbind(data.frame(tax_table(tmpc1), rowSums(otu_table(tmpc1))), "HI_worm_1428")
colnames(tmptaxc)[8] <- "ReadSums"
colnames(tmptaxc)[9] <- "MTD"

tax1428 <- rbind(tmptaxa, tmptaxb, tmptaxc)
#


#   ##### For D28 to D41 #####

lotreatearly <- ps28lotreat
loctrlearly <- ps28loctrl
hitreatearly <- ps28hitreat

lotreattax <- as.list(rownames(otu_table(lotreatearly)))
loctrltax <- as.list(rownames(otu_table(loctrlearly)))
hitreattax <- as.list(rownames(otu_table(hitreatearly)))

lotreatlate <- ps41lotreat
loctrllate <- ps41loctrl
hitreatlate <- ps41hitreat

length(loctrltax) #
tmpa <- subset(otu_table(lotreatlate), !rownames(otu_table(lotreatlate)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatlate), metadata=sample_data(lotreatlate))
#28 taxa remaining
richa <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "LO_worm_2841",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpa1)), 
                     "Family" = data.frame(tax_table(tmpa1))$Family)
richa1 <- merge(richa, tmpdat, by="ASV", all=FALSE)

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrllate), !rownames(otu_table(loctrllate)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrllate), metadata=sample_data(loctrllate))
#17 taxa remaining

tmpc <- subset(otu_table(hitreatlate), !rownames(otu_table(hitreatlate)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatlate), metadata=sample_data(hitreatlate))
#171 taxa remaining
richc <- ResistDiffAbundBC[ResistDiffAbundBC$MTD == "HI_worm_2841",]
tmpdat <- data.frame(ASV = rownames(tax_table(tmpc1)), 
                     "Family" = data.frame(tax_table(tmpc1))$Family)
richc1 <- merge(richc, tmpdat, by="ASV", all=FALSE)

#Now let's bring everything together into one large object. 
ps2841invaders <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tree = phy_tree(psexp3000))
phyobj <- ps2841invaders
sample_data(phyobj)$comparison <- "1428"
sample_data(phyobj)$readsInvaders <- sample_sums(phyobj)
tmp <- data.frame(estimate_richness(phyobj, measures="Observed"))
sample_data(phyobj)$richnessInvaders <- tmp$Observed
sample_data(phyobj)$sample <- rownames(sample_data(phyobj))
ps2841invaders <- phyobj
#View(sample_data(ps2841invaders))

#Now extract the taxa from the invaders and quantify, so that I can merge with the 
tmptaxa <- cbind(data.frame(tax_table(tmpa1), rowSums(otu_table(tmpa1))), "LO_worm_2841")
colnames(tmptaxa)[8] <- "ReadSums"
colnames(tmptaxa)[9] <- "MTD"
tmptaxb <- cbind(data.frame(tax_table(tmpb1), rowSums(otu_table(tmpb1))), "LO_sham_2841")
colnames(tmptaxb)[8] <- "ReadSums"
colnames(tmptaxb)[9] <- "MTD"
tmptaxc <- cbind(data.frame(tax_table(tmpc1), rowSums(otu_table(tmpc1))), "HI_worm_2841")
colnames(tmptaxc)[8] <- "ReadSums"
colnames(tmptaxc)[9] <- "MTD"

tax2841 <- rbind(tmptaxa, tmptaxb, tmptaxc)
#
#
#     ##### Combine data into a single phyloseq object #####
otu_table <- phyloseq::otu_table
tax_table <- phyloseq::tax_table

#Keeps erroring, I think due to different OTU and Tax tables. Would need to pull 
# the phyloseq objects apart and rebuild. 
#
psallinvaders <- merge_phyloseq(otu_table(pspre1invaders), tax_table(pspre1invaders), sample_data(pspre1invaders), 
               otu_table(ps110invaders), tax_table(ps110invaders), sample_data(ps110invaders), 
               otu_table(ps1014invaders), tax_table(ps1014invaders), sample_data(ps1014invaders), 
               otu_table(ps1428invaders), tax_table(ps1428invaders), sample_data(ps1428invaders), 
               otu_table(ps2841invaders), tax_table(ps2841invaders), sample_data(ps2841invaders)) 
               
#
#     ##### Combine data into a single df to assess whether things differed between groups #####
dftmp <- data.frame(rbind(sample_data(pspre1invaders), sample_data(ps110invaders), sample_data(ps1014invaders), 
               sample_data(ps1428invaders), sample_data(ps2841invaders)))
dftmp <- data.frame(sample_data(psallinvaders))
#Pull only the columns I want for plotting and future analysis. 
dfinvaderdat <- dftmp[,c("sample", "richnessInvaders", "readsInvaders", "mouse", "MT", "MTD", "microb", "day", "day_num", "comparison", "treatment")]
#


#     ##### Create stacked bars for each of the average invader taxa by group #####
#       ##### Family colors #####
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

#       ##### Create a stacked bar plot to check out family=level invaders #####
#View(sample_data(psallinvaders))
psfam <- aggregate_taxa(psallinvaders, level="Family")
pseq1 <- aggregate_rare(psallinvaders, level = "Family", detection = 0.01, prevalence = 0.05)

nrow(unique(tax_table(pseq1)[, "Family"])) #15 families when aggregated rare for psallinvaders

library(ggbreak)

ppinvaders <- plot_composition(pseq1, 
                               taxonomic.level="Family", 
                               otu.sort= "abundance", 
                               #transform = "compositional",
                               average_by = "MTD", 
                               plot_type='barplot', 
                               group_by="MT") +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=famcolsall) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=1),
        plot.title=element_text(size=13, face="bold"), 
        axis.title=element_text(size=13, face="bold"), 
        axis.text=element_text(size=12), 
        legend.title=element_text(size=11, face="bold"), 
        legend.text=element_text(size=10)) + 
  scale_y_break(c(150, 375), ticklabels = c(400), scales="fixed", expand=TRUE)+
  guides(fill=guide_legend(title="Family")) +
  labs(x = "Group", y = "Read Abundance (out of 3000)", title="Read counts Invading bacteria") + 
  ggpubr::rotate_x_text(40) 
ppinvaders
#
#
#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/AllInvadingASVs_StackedBar_byGroup.pdf", 
#       ppinvaders, units="in", height=6, width=7)
#
#   ##### Plot average richness of invaders over time by group #####
dfhitreat <- dfinvaderdat[dfinvaderdat$microb=="HI" & dfinvaderdat$treatment=="worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(richnessInvaders)
dfhictrl <- dfinvaderdat[dfinvaderdat$microb=="HI" & dfinvaderdat$treatment=="control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(richnessInvaders)
dflotreat <- dfinvaderdat[dfinvaderdat$microb=="LO" & dfinvaderdat$treatment=="worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(richnessInvaders)
dfloctrl <- dfinvaderdat[dfinvaderdat$microb=="LO" & dfinvaderdat$treatment=="control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(richnessInvaders)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessInvaders = mean(richnessInvaders, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessInvaders = mean(richnessInvaders, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessInvaders = mean(richnessInvaders, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessInvaders = mean(richnessInvaders, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessInvaders, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessInvaders, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessInvaders, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessInvaders, na.rm=TRUE))

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

dfhitreat.mean$SEH <- dfhitreat.mean$meanrichnessInvaders + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanrichnessInvaders + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanrichnessInvaders + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanrichnessInvaders + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanrichnessInvaders - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanrichnessInvaders - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanrichnessInvaders - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanrichnessInvaders - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as ASVrichnessInvaders-specific dfs. 
dfhitreat.mean.ASVrichnessInvaders <- dfhitreat.mean
##  day_num meanrichnessInvaders    SD     N   SEH    SEL
##1       1                 4.28  3.16    18  5.02  3.53 
##2      10                 4.44  3.24    18  5.21  3.68 
##3      14                 6.07  3.63    15  7.01  5.13 
##4      28                 1.6   2.51     5  2.72  0.478
##5      41                 40.5  21.6      8 48.1  32.9
dfhictrl.mean.ASVrichnessInvaders <- dfhictrl.mean 
##1       1                 16.5 11.2      4  22.1  10.9
##2      10                 39.5 13.5      4  46.2  32.8
##3      14                 15.2  9.54     4  20.0  10.5
dflotreat.mean.ASVrichnessInvaders <- dflotreat.mean
##       1                 1.09   2.07     32  1.46   0.728
##2      10                28.1   11.4     32 30.1   26.1  
##3      14                0.531  0.915    32  0.693  0.369
##4      28                0.867  0.834    15  1.08   0.651
##5      41                3.70   1.79     23  4.07   3.32 
dfloctrl.mean.ASVrichnessInvaders <- dfloctrl.mean 
##1       1                 3.33  2.58     6  4.39  2.28
##2      10                 4.17  3.06     6  5.42  2.92
##3      14                 2.67  2.73     6  3.78  1.55
##4      28                 3.67  3.98     6  5.29  2.04
##5      41                 4.33  2.34     6  5.29  3.38

#         ##### Plot as average over time with error bars #####
ASVrichnessInvaders <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfhitreat.mean.ASVrichnessInvaders$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfhictrl.mean.ASVrichnessInvaders$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.ASVrichnessInvaders$SEL, ymax=dfhitreat.mean.ASVrichnessInvaders$SEH, y=dfhitreat.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfhitreat.mean.ASVrichnessInvaders$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.ASVrichnessInvaders$SEL, ymax=dfhictrl.mean.ASVrichnessInvaders$SEH, y=dfhictrl.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfhictrl.mean.ASVrichnessInvaders$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfhitreat.mean.ASVrichnessInvaders$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfhictrl.mean.ASVrichnessInvaders$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dflotreat.mean.ASVrichnessInvaders$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfloctrl.mean.ASVrichnessInvaders$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.ASVrichnessInvaders$SEL, ymax=dflotreat.mean.ASVrichnessInvaders$SEH, y=dflotreat.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dflotreat.mean.ASVrichnessInvaders$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.ASVrichnessInvaders$SEL, ymax=dfloctrl.mean.ASVrichnessInvaders$SEH, y=dfloctrl.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfloctrl.mean.ASVrichnessInvaders$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dflotreat.mean.ASVrichnessInvaders$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.ASVrichnessInvaders$meanrichnessInvaders, x=dfloctrl.mean.ASVrichnessInvaders$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("ASV richness") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average ASV richness Invaders") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

ASVrichnessInvaders

#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_AverageASVrichnessInvaders_Resistance_lineplot_woD31.pdf", ASVrichnessInvaders, units="in", height=5, width=6)
#




#
#   ##### Plot average readsInvaders of invaders over time by group #####

dfhitreat <- dfinvaderdat[dfinvaderdat$microb=="HI" & dfinvaderdat$treatment=="worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(readsInvaders)
dfhictrl <- dfinvaderdat[dfinvaderdat$microb=="HI" & dfinvaderdat$treatment=="control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(readsInvaders)
dflotreat <- dfinvaderdat[dfinvaderdat$microb=="LO" & dfinvaderdat$treatment=="worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(readsInvaders)
dfloctrl <- dfinvaderdat[dfinvaderdat$microb=="LO" & dfinvaderdat$treatment=="control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(readsInvaders)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsInvaders = mean(readsInvaders, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsInvaders = mean(readsInvaders, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsInvaders = mean(readsInvaders, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsInvaders = mean(readsInvaders, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsInvaders, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsInvaders, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsInvaders, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsInvaders, na.rm=TRUE))

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

dfhitreat.mean$SEH <- dfhitreat.mean$meanreadsInvaders + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanreadsInvaders + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanreadsInvaders + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanreadsInvaders + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanreadsInvaders - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanreadsInvaders - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanreadsInvaders - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanreadsInvaders - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as ASVreadsInvaders-specific dfs. 
dfhitreat.mean.ASVreadsInvaders <- dfhitreat.mean
## day_num meanrichness    SD     N   SEH    SEL
##1       1         4.28  3.16    18  5.02  3.53 
##2      10         4.44  3.24    18  5.21  3.68 
##3      14         6.07  3.63    15  7.01  5.13 
##4      28         1.6   2.51     5  2.72  0.478
##5      41        40.5  21.6      8 48.1  32.9  
dfhictrl.mean.ASVreadsInvaders <- dfhictrl.mean 
##1       1      32.2  30.0     4  47.3  17.2
##2      10      73.2  34.4     4  90.4  56.1
##3      14      36    20.8     4  46.4  25.6
dflotreat.mean.ASVreadsInvaders <- dflotreat.mean
##1       1      9.88  36.5     32  16.3    3.43 
##2      10      439.   289.      32 490.   388.   
##3      14      1.12   2.04    32   1.49   0.764
##4      28     10     18.7     15  14.8    5.16 
##5      41     19.9   18.7     23  23.8   16.0  
dfloctrl.mean.ASVreadsInvaders <- dfloctrl.mean 
##1       1       7.5  9.14     6  11.2  3.77
##2      10      35.7 27.9      6  47.1 24.3 
##3      14       8.5 10.7      6  12.9  4.14
##4      28       9   12.2      6  14.0  4.01
##5      41      13.8  7.65     6  17.0 10.7


#         ##### Plot as average over time with error bars #####
ASVreadsInvaders <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.ASVreadsInvaders$meanreadsInvaders, x=dfhitreat.mean.ASVreadsInvaders$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.ASVreadsInvaders$meanreadsInvaders, x=dfhictrl.mean.ASVreadsInvaders$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.ASVreadsInvaders$SEL, ymax=dfhitreat.mean.ASVreadsInvaders$SEH, y=dfhitreat.mean.ASVreadsInvaders$meanreadsInvaders, x=dfhitreat.mean.ASVreadsInvaders$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.ASVreadsInvaders$SEL, ymax=dfhictrl.mean.ASVreadsInvaders$SEH, y=dfhictrl.mean.ASVreadsInvaders$meanreadsInvaders, x=dfhictrl.mean.ASVreadsInvaders$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.ASVreadsInvaders$meanreadsInvaders, x=dfhitreat.mean.ASVreadsInvaders$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.ASVreadsInvaders$meanreadsInvaders, x=dfhictrl.mean.ASVreadsInvaders$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.ASVreadsInvaders$meanreadsInvaders, x=dflotreat.mean.ASVreadsInvaders$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.ASVreadsInvaders$meanreadsInvaders, x=dfloctrl.mean.ASVreadsInvaders$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.ASVreadsInvaders$SEL, ymax=dflotreat.mean.ASVreadsInvaders$SEH, y=dflotreat.mean.ASVreadsInvaders$meanreadsInvaders, x=dflotreat.mean.ASVreadsInvaders$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.ASVreadsInvaders$SEL, ymax=dfloctrl.mean.ASVreadsInvaders$SEH, y=dfloctrl.mean.ASVreadsInvaders$meanreadsInvaders, x=dfloctrl.mean.ASVreadsInvaders$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.ASVreadsInvaders$meanreadsInvaders, x=dflotreat.mean.ASVreadsInvaders$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.ASVreadsInvaders$meanreadsInvaders, x=dfloctrl.mean.ASVreadsInvaders$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("Number readsInvaders (out of 3000)") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average readsInvaders") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

ASVreadsInvaders

#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_AverageReadsInvaders_Resistance_lineplot_woD31.pdf", ASVreadsInvadersInvaders, units="in", height=5, width=6)
#

#

# ##### Create a boxplot for the invaders over time instead of line plot #####

View(dfinvaderdat)
dfinvaderdat$MT <- factor(dfinvaderdat$MT, levels=c("LO_worm", "LO_control", "HI_worm", "HI_control"))
dfinvaderdat$MTD <- factor(dfinvaderdat$MTD, 
                           levels=c('LO_worm_D1', 'LO_worm_D10' , 'LO_worm_D14', 'LO_worm_D28' , 'LO_worm_D41' ,
                                    'LO_control_D1', 'LO_control_D10','LO_control_D14', 'LO_control_D28', 'LO_control_D41' ,
                                    'HI_worm_D1' ,  'HI_worm_D10','HI_worm_D14',  'HI_worm_D28','HI_worm_D41',
                                    'HI_control_D1' ,'HI_control_D10', 'HI_control_D14'))
allcols2 <- c('#c275eb' , '#A11FE5', '#7E05BE' , '#5C008D', '#300049',
              '#c275eb' , '#A11FE5', '#7E05BE' , '#5C008D', '#300049',   
              "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" ,"#644C00",
              "#f2d06d" , "#DCAE1C" , "#BC8F00" )


allfill1 <- c('#c275eb' , '#A11FE5' , '#7E05BE' , '#5C008D', '#300049',
              "#FFFFFF",  "#FFFFFF",   "#FFFFFF", "#FFFFFF", "#FFFFFF",
              "#f2d06d" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00" ,
              "#FFFFFF",  "#FFFFFF",  "#FFFFFF")


invaderbox <- ggplot(dfinvaderdat, aes(x=day_num, y=richnessInvaders)) +
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  geom_boxplot(aes(color=MTD), width=2, position = position_dodge(width=2), outlier.shape = NA) +
  geom_point(aes(shape=MT, color=MTD, fill=MTD), position = position_jitterdodge(jitter.width = 0.5, dodge.width = 2), alpha=0.9, size=3) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_color_manual(values=c(allcols2)) +
  scale_fill_manual(values=c(allfill1)) +
  theme_bw() +
  ylab("Number of Invaders") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average Invader richness") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
invaderbox

ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/InvadersBoxplotAllDays_woD31.pdf", 
       invaderbox, units="in", height=5, width=6)
  

#   ##### Write out the concatenated OTU and taxa tables for these "invader" phyloseq objects. ##### 
OTU <- data.frame(otu_table(pspre1invaders))
OTU$sum <- rowSums(OTU)
pspre1invader_otutax <- cbind(tax_table(pspre1invaders), OTU, "ASV"=rownames(OTU))
pspre1invader_otutax$comparison <- "pre1"
View(pspre1invader_otutax)

OTU <- data.frame(otu_table(ps110invaders))
OTU$sum <- rowSums(OTU)
ps110invader_otutax <- cbind(tax_table(ps110invaders), OTU, "ASV"=rownames(OTU))
ps110invader_otutax$comparison <- "110"
View(ps110invader_otutax)

OTU <- data.frame(otu_table(ps1014invaders))
OTU$sum <- rowSums(OTU)
ps1014invader_otutax <- cbind(tax_table(ps1014invaders), OTU, "ASV"=rownames(OTU))
ps1014invader_otutax$comparison <- "1014"
View(ps1014invader_otutax)

OTU <- data.frame(otu_table(ps1428invaders))
OTU$sum <- rowSums(OTU)
ps1428invader_otutax <- cbind(tax_table(ps1428invaders), OTU, "ASV"=rownames(OTU))
ps1428invader_otutax$comparison <- "1428"
View(ps1428invader_otutax)

OTU <- data.frame(otu_table(ps2841invaders))
OTU$sum <- rowSums(OTU)
ps2841invader_otutax <- cbind(tax_table(ps2841invaders), OTU, "ASV"=rownames(OTU))
ps2841invader_otutax$comparison <- "2841"
View(ps2841invader_otutax)


write.csv(pre1invader_otutax, "~/Documents/1U4U_16S/GG2_reclassification/OTU_Tax_tables/pspre1invader_otutax.csv")
write.csv(ps110invader_otutax, "~/Documents/1U4U_16S/GG2_reclassification/OTU_Tax_tables/ps110invader_otutax.csv")
write.csv(ps1014invader_otutax, "~/Documents/1U4U_16S/GG2_reclassification/OTU_Tax_tables/ps1014invader_otutax.csv")
write.csv(ps1428invader_otutax, "~/Documents/1U4U_16S/GG2_reclassification/OTU_Tax_tables/ps1428invader_otutax.csv")
write.csv(ps2841invader_otutax, "~/Documents/1U4U_16S/GG2_reclassification/OTU_Tax_tables/ps2841invader_otutax.csv")

#

##### FIND WHETHER THE INVADERS WERE ALSO DIFFERENTIALLY ABUNDANT #####
head(ResistDiffAbundBC)
nrow(ResistDiffAbundBC) #This is the DF we want to use after running ANCOM-BC. 

#Binding all the rows of the df's that I made from the invading taxa tables. 
invaderTax <- rbind(taxpre1, tax110, tax1014, tax1428, tax2841)
#Check it out
head(invaderTax)
#Create the columns that will help relate the diff abund and invader dfs to each other. 
invaderTax$shortASV <- strtrim(rownames(invaderTax), 22)
invaderTax$mergecol <- paste(invaderTax$MTD, invaderTax$shortASV, sep="_")

# Create a new column 'invading' based on whether the same mergecol is in both df's
ResistDiffAbundBC$invading <- ifelse(ResistDiffAbundBC$mergecol %in% invaderTax$mergecol, "yes", "no")
#
#



#
##### FINDING ASVs LOST BY EACH GROUP #####
# Going to use the same approach as invaders, but going to subtract the late from the early ######
#   ##### For pretreatment to D1 #####

lotreatearly <- psprelotreat
loctrlearly <- pspreloctrl
hitreatearly <- psprehitreat
hictrlearly <- psprehictrl

lotreatlate <- ps1lotreat
loctrllate <- ps1loctrl
hitreatlate <- ps1hitreat
hictrllate <- ps1hictrl

lotreattax <- as.list(rownames(otu_table(lotreatlate)))
loctrltax <- as.list(rownames(otu_table(loctrllate)))
hitreattax <- as.list(rownames(otu_table(hitreatlate)))
hictrltax <- as.list(rownames(otu_table(hictrllate)))

length(lotreattax) #110 ASVs
tmpa <- subset(otu_table(lotreatearly), !rownames(otu_table(lotreatearly)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatearly), metadata=sample_data(lotreatearly))
#19 taxa 

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrlearly), !rownames(otu_table(loctrlearly)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrlearly), metadata=sample_data(loctrlearly))
#8 taxa remaining

tmpc <- subset(otu_table(hitreatearly), !rownames(otu_table(hitreatearly)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatearly), metadata=sample_data(hitreatearly))
#42 taxa lost

tmpd <- subset(otu_table(hictrlearly), !rownames(otu_table(hictrlearly)) %in% hictrltax)
tmpd1 <- merge_phyloseq(tmpd, taxonomy=tax_table(hictrlearly), metadata=sample_data(hictrlearly))
#65 taxa remaining

#Now let's bring everything together into one large object. 
pspre1lost <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tmpd1, tree = phy_tree(psexp3000))
sample_data(pspre1lost)$readsLost <- sample_sums(pspre1lost)
sample_data(pspre1lost)$comparison <- "pre1"
#123 total ASVs in the object. 
#

#   ##### For D1 to D10 #####

lotreatearly <- ps1lotreat
loctrlearly <- ps1loctrl
hitreatearly <- ps1hitreat
hictrlearly <- ps1hictrl

lotreatlate <- ps10lotreat
loctrllate <- ps10loctrl
hitreatlate <- ps10hitreat
hictrllate <- ps10hictrl

lotreattax <- as.list(rownames(otu_table(lotreatlate)))
loctrltax <- as.list(rownames(otu_table(loctrllate)))
hitreattax <- as.list(rownames(otu_table(hitreatlate)))
hictrltax <- as.list(rownames(otu_table(hictrllate)))


length(lotreattax) #89 ASVs
tmpa <- subset(otu_table(lotreatearly), !rownames(otu_table(lotreatearly)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatearly), metadata=sample_data(lotreatearly))
#11 ASVs

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrlearly), !rownames(otu_table(loctrlearly)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrlearly), metadata=sample_data(loctrlearly))
#10 ASVs

tmpc <- subset(otu_table(hitreatearly), !rownames(otu_table(hitreatearly)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatearly), metadata=sample_data(hitreatearly))
#50 ASVs

tmpd <- subset(otu_table(hictrlearly), !rownames(otu_table(hictrlearly)) %in% hictrltax)
tmpd1 <- merge_phyloseq(tmpd, taxonomy=tax_table(hictrlearly), metadata=sample_data(hictrlearly))
#34 ASVs

#Now let's bring everything together into one large object. 
ps110lost <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tmpd1, tree = phy_tree(psexp3000no31))
sample_data(ps110lost)$readsLost <- sample_sums(ps110lost)
sample_data(ps110lost)$comparison <- "110"
#95 total ASVs

#   ##### For D10 to D14 #####

lotreatearly <- ps10lotreat
loctrlearly <- ps10loctrl
hitreatearly <- ps10hitreat
hictrlearly <- ps10hictrl

lotreatlate <- ps14lotreat
loctrllate <- ps14loctrl
hitreatlate <- ps14hitreat
hictrllate <- ps14hictrl

lotreattax <- as.list(rownames(otu_table(lotreatlate)))
loctrltax <- as.list(rownames(otu_table(loctrllate)))
hitreattax <- as.list(rownames(otu_table(hitreatlate)))
hictrltax <- as.list(rownames(otu_table(hictrllate)))


length(loctrltax) #105 ASVs
tmpa <- subset(otu_table(lotreatearly), !rownames(otu_table(lotreatearly)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatearly), metadata=sample_data(lotreatearly))
tmpa1
#33taxa remaining

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrlearly), !rownames(otu_table(loctrlearly)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrlearly), metadata=sample_data(loctrlearly))
#14 taxa remaining

tmpc <- subset(otu_table(hitreatearly), !rownames(otu_table(hitreatearly)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatearly), metadata=sample_data(hitreatearly))
#48 taxa remaining

tmpd <- subset(otu_table(hictrlearly), !rownames(otu_table(hictrlearly)) %in% hictrltax)
tmpd1 <- merge_phyloseq(tmpd, taxonomy=tax_table(hictrlearly), metadata=sample_data(hictrlearly))
#70 taxa remaining

#Now let's bring everything together into one large object. 
ps1014lost <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tmpd1, tree = phy_tree(psexp3000))
sample_data(ps1014lost)$readsLost <- sample_sums(ps1014lost)
sample_data(ps1014lost)$comparison <- "1014"
#153 ASVs 

#   ##### For D14 to D28 #####

lotreatearly <- ps14lotreat
loctrlearly <- ps14loctrl
hitreatearly <- ps14hitreat
hictrlearly <- ps14hictrl

lotreatlate <- ps28lotreat
loctrllate <- ps28loctrl
hitreatlate <- ps28hitreat

lotreattax <- as.list(rownames(otu_table(lotreatlate)))
loctrltax <- as.list(rownames(otu_table(loctrllate)))
hitreattax <- as.list(rownames(otu_table(hitreatlate)))

length(lotreattax) #
tmpa <- subset(otu_table(lotreatearly), !rownames(otu_table(lotreatearly)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatearly), metadata=sample_data(lotreatearly))
#44 taxa remaining

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrlearly), !rownames(otu_table(loctrlearly)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrlearly), metadata=sample_data(loctrlearly))
#14 taxa remaining

tmpc <- subset(otu_table(hitreatearly), !rownames(otu_table(hitreatearly)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatearly), metadata=sample_data(hitreatearly))
#240 taxa remaining

#Now let's bring everything together into one large object. 
ps1428lost <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tree = phy_tree(psexp3000))
sample_data(ps1428lost)$readsLost <- sample_sums(ps1428lost)
sample_data(ps1428lost)$comparison <- "1428"
#281 ASVs in whole dataset

#   ##### For D28 to D41 #####

lotreatearly <- ps28lotreat
loctrlearly <- ps28loctrl
hitreatearly <- ps28hitreat

lotreatlate <- ps41lotreat
loctrllate <- ps41loctrl
hitreatlate <- ps41hitreat
hictrllate <- ps41hictrl

lotreattax <- as.list(rownames(otu_table(lotreatlate)))
loctrltax <- as.list(rownames(otu_table(loctrllate)))
hitreattax <- as.list(rownames(otu_table(hitreatlate)))
hictrltax <- as.list(rownames(otu_table(hictrllate)))

length(loctrltax) #
tmpa <- subset(otu_table(lotreatearly), !rownames(otu_table(lotreatearly)) %in% lotreattax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(lotreatearly), metadata=sample_data(lotreatearly))
#12 taxa remaining

#Do this for the other phyloseq objects. 
tmpb <- subset(otu_table(loctrlearly), !rownames(otu_table(loctrlearly)) %in% loctrltax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(loctrlearly), metadata=sample_data(loctrlearly))
#15 taxa remaining

tmpc <- subset(otu_table(hitreatearly), !rownames(otu_table(hitreatearly)) %in% hitreattax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(hitreatearly), metadata=sample_data(hitreatearly))
#17 taxa remaining

#Now let's bring everything together into one large object. 
ps2841lost <- merge_phyloseq(tmpa1, tmpb1, tmpc1, tree = phy_tree(psexp3000))
sample_data(ps2841lost)$readsLost <- sample_sums(ps2841lost)
sample_data(ps2841lost)$comparison <- "2841"
#41 ASVs total




#

#   ##### Create concatenated phloseq objects and sample data df for lost ASVs ######
#into a phyloseq object
psalllost <- merge_phyloseq(otu_table(pspre1lost), tax_table(pspre1lost), sample_data(pspre1lost), 
                                otu_table(ps110lost), tax_table(ps110lost), sample_data(ps110lost), 
                                otu_table(ps1014lost), tax_table(ps1014lost), sample_data(ps1014lost), 
                                otu_table(ps1428lost), tax_table(ps1428lost), sample_data(ps1428lost), 
                                otu_table(ps2841lost), tax_table(ps2841lost), sample_data(ps2841lost)) 

tmp <- estimate_richness(psalllost, measures="Observed")
sample_data(psalllost)$richnessLost <- tmp$Observed
sample_data(psalllost)$mouse_comp <- paste(sample_data(psalllost)$mouse, sample_data(psalllost)$comparison, sep="_")


#Now pull it into a dataframe
dflostdatall <- data.frame(sample_data(psalllost))
#Pull only the columns I want for plotting and future analysis. Going to merge this with dfinvaderdat. 
dflostdat <- dflostdatall[,c("richnessLost", "readsLost", "mouse_comp")]

#
#

##### FINDING THE INVADING ASVs THAT COULD HAVE BEEN SOURCED FROM OPPOSITE MICROBIOME #####
#Looking for the intersection of the invading ASVs from the current timepoint of 1 microbiome type
# and the ASVs from the other microbiome type at all the timepoints before that one. 


#   ##### Create phyloseq objects that have accruing taxa from all timepoints instead of just 1 or two timepoints. ##### 
#Pull them together from their pieces since they have different OTUs in them. 
#Pretreatment only is already made. 

#Make those with pre and 1 dpi. 
pshipre1 <- merge_phyloseq(tax_table(psprehitreat),tax_table(ps1hitreat), 
                             otu_table(psprehitreat), otu_table(ps1hitreat), 
                             sample_data(psprehitreat), sample_data(ps1hitreat), 
                             tax_table(psprehictrl),tax_table(ps1hictrl),
                             otu_table(psprehictrl), otu_table(ps1hictrl), 
                             sample_data(psprehictrl), sample_data(ps1hictrl))

pslopre1 <- merge_phyloseq(tax_table(psprelotreat),tax_table(ps1lotreat), 
                             otu_table(psprelotreat), otu_table(ps1lotreat), 
                             sample_data(psprelotreat), sample_data(ps1lotreat), 
                             tax_table(pspreloctrl),tax_table(ps1loctrl), 
                             otu_table(pspreloctrl), otu_table(ps1loctrl), 
                             sample_data(pspreloctrl), sample_data(ps1loctrl))

#Make those with pre, 1, 10 dpi. 
pshipre110 <- merge_phyloseq(tax_table(psprehitreat),tax_table(ps1hitreat), tax_table(ps10hitreat),
                                 otu_table(psprehitreat), otu_table(ps1hitreat), tax_table(ps10hitreat), 
                                 sample_data(psprehitreat), sample_data(ps1hitreat), sample_data(ps10hitreat), 
                                 tax_table(psprehictrl),tax_table(ps1hictrl), tax_table(ps10hictrl),
                                 otu_table(psprehictrl), otu_table(ps1hictrl), tax_table(ps10hictrl), 
                                 sample_data(psprehictrl), sample_data(ps1hictrl), sample_data(ps10hictrl))

pslopre110 <- merge_phyloseq(tax_table(psprelotreat),tax_table(ps1lotreat), tax_table(ps10lotreat),
                                 otu_table(psprelotreat), otu_table(ps1lotreat), tax_table(ps10lotreat), 
                                 sample_data(psprelotreat), sample_data(ps1lotreat), sample_data(ps10lotreat), 
                                 tax_table(pspreloctrl),tax_table(ps1loctrl), tax_table(ps10loctrl),
                                 otu_table(pspreloctrl), otu_table(ps1loctrl), tax_table(ps10loctrl), 
                                 sample_data(pspreloctrl), sample_data(ps1loctrl), sample_data(ps10loctrl))
#
#Make those with pre, 1, 10, 14 dpi. 
pshipre11014 <- merge_phyloseq(tax_table(psprehitreat),tax_table(ps1hitreat), tax_table(ps10hitreat),
                                 otu_table(psprehitreat), otu_table(ps1hitreat), tax_table(ps10hitreat), 
                                 sample_data(psprehitreat), sample_data(ps1hitreat), sample_data(ps10hitreat), 
                                 tax_table(ps14hitreat),otu_table(ps14hitreat),  sample_data(ps14hitreat),
                                 tax_table(psprehictrl),tax_table(ps1hictrl), tax_table(ps10hictrl),
                                 otu_table(psprehictrl), otu_table(ps1hictrl), tax_table(ps10hictrl), 
                                 sample_data(psprehictrl), sample_data(ps1hictrl), sample_data(ps10hictrl), 
                                 tax_table(ps14hictrl),otu_table(ps14hictrl),  sample_data(ps14hictrl))

pslopre11014 <- merge_phyloseq(tax_table(psprelotreat),tax_table(ps1lotreat), tax_table(ps10lotreat),
                                 otu_table(psprelotreat), otu_table(ps1lotreat), tax_table(ps10lotreat), 
                                 sample_data(psprelotreat), sample_data(ps1lotreat), sample_data(ps10lotreat), 
                                 tax_table(ps14lotreat),otu_table(ps14lotreat),  sample_data(ps14lotreat),
                                 tax_table(pspreloctrl),tax_table(ps1loctrl), tax_table(ps10loctrl),
                                 otu_table(pspreloctrl), otu_table(ps1loctrl), tax_table(ps10loctrl), 
                                 sample_data(pspreloctrl), sample_data(ps1loctrl), sample_data(ps10loctrl), 
                                 tax_table(ps14loctrl),otu_table(ps14loctrl),  sample_data(ps14loctrl))
#
#All but 41 dpi
pshipre1101428 <- merge_phyloseq(tax_table(psprehitreat),tax_table(ps1hitreat), tax_table(ps10hitreat),
                             otu_table(psprehitreat), otu_table(ps1hitreat), tax_table(ps10hitreat), 
                             sample_data(psprehitreat), sample_data(ps1hitreat), sample_data(ps10hitreat), 
                             tax_table(ps14hitreat),otu_table(ps14hitreat),  sample_data(ps14hitreat),
                             tax_table(ps28hitreat),otu_table(ps28hitreat),  sample_data(ps28hitreat),
                             tax_table(psprehictrl),tax_table(ps1hictrl), tax_table(ps10hictrl),
                             otu_table(psprehictrl), otu_table(ps1hictrl), tax_table(ps10hictrl), 
                             sample_data(psprehictrl), sample_data(ps1hictrl), sample_data(ps10hictrl), 
                             tax_table(ps14hictrl),otu_table(ps14hictrl),  sample_data(ps14hictrl))

pslopre1101428 <- merge_phyloseq(tax_table(psprelotreat),tax_table(ps1lotreat), tax_table(ps10lotreat),
                                 otu_table(psprelotreat), otu_table(ps1lotreat), tax_table(ps10lotreat), 
                                 sample_data(psprelotreat), sample_data(ps1lotreat), sample_data(ps10lotreat), 
                                 tax_table(ps14lotreat),otu_table(ps14lotreat),  sample_data(ps14lotreat),
                                 tax_table(ps28lotreat),otu_table(ps28lotreat),  sample_data(ps28lotreat),
                                 tax_table(pspreloctrl),tax_table(ps1loctrl), tax_table(ps10loctrl),
                                 otu_table(pspreloctrl), otu_table(ps1loctrl), tax_table(ps10loctrl), 
                                 sample_data(pspreloctrl), sample_data(ps1loctrl), sample_data(ps10loctrl), 
                                 tax_table(ps14loctrl),otu_table(ps14loctrl),  sample_data(ps14loctrl),
                                 tax_table(ps28loctrl),otu_table(ps28loctrl),  sample_data(ps28loctrl))
#


#   ##### First off, just find the number of shared ASVs between LO and HI at pretreatment.  #####
lotax <- as.list(rownames(otu_table(psprelo)))
hitax <- as.list(rownames(otu_table(psprehi)))

#Number of HI taxa in LO at pre
HIinLO <- subset(otu_table(psprelo), rownames(otu_table(psprelo)) %in% hitax)
pspreHIinLO <- merge_phyloseq(HIinLO, taxonomy=tax_table(psprelo), metadata=sample_data(psprelo))

dfHIinLO <- data.frame(sample <- rownames(sample_data(pspreHIinLO)), 
                       PreOverlapReads <- sample_sums(pspreHIinLO), 
                       PreOverlapRichness  <- estimate_richness(pspreHIinLO, measures="Observed"),
                       treatment = sample_data(pspreHIinLO)$treatment)
#Rename the columns. 
colnames(dfHIinLO)[1] <- "sample"
colnames(dfHIinLO)[2] <- "PreOverlapReads"
colnames(dfHIinLO)[3] <- "PreOverlapRichness"
dfHIinLO$microb <- "LO"


#Number of LO taxa in HI at pre
LOinHI <- subset(otu_table(psprehi), rownames(otu_table(psprehi)) %in% lotax)
pspreLOinHI <- merge_phyloseq(LOinHI, taxonomy=tax_table(psprehi), metadata=sample_data(psprehi))

dfLOinHI <- data.frame(sample <- rownames(sample_data(pspreLOinHI)), 
                       PreOverlapReads <- sample_sums(pspreLOinHI), 
                       PreOverlapRichness  <- estimate_richness(pspreLOinHI, measures="Observed"),
                       treatment = sample_data(pspreLOinHI)$treatment)
colnames(dfLOinHI)[1] <- "sample"
colnames(dfLOinHI)[2] <- "PreOverlapReads"
colnames(dfLOinHI)[3] <- "PreOverlapRichness"
dfLOinHI$microb <- "HI"
head(dfLOinHI)
#
#Combine these pretreatment overlap read dfs.
dfpreOverlap <- rbind(dfLOinHI, dfHIinLO)
dfpreOverlap$MT <- paste(dfpreOverlap$microb, dfpreOverlap$treatment, sep="_")
#
#   ##### Now calculate the number of ASVs that could have come from the other microbiome  #####
#         type previous to the target timepoint vs. those that are not. 
#      ##### pretreatment to D1. #####
#starting with objects that only have pretreatment samples. 
psprelo
psprehi

#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(pspre1invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#93 taxa here in HI out of the total 122
psloinvaders <- subset_samples(pspre1invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#30 taxa here in LO out of the total 122
#Meaning they overlap in 1 ASV. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(psprelo)))
hitax <- as.list(rownames(otu_table(psprehi)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% lotax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#13 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% hitax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
#22 taxa remaining 
#
pspre1sharedtax <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                  otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(pspre1sharedtax)$comparison <- "pre1"
#
                                  
#
#      ##### D1 to D10 #####

#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps110invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#148 taxa here in HI out of the total 261
psloinvaders <- subset_samples(ps110invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#137 taxa here in LO out of the total 261. 
#Meaning they overlap in 24 ASVs. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre1)))
hitax <- as.list(rownames(otu_table(pshipre1)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% lotax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#21 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% hitax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
#84  taxa remaining 
#
ps110sharedtax <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                 otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps110sharedtax)$comparison <- "110"
#
#
#      ##### D10 to D14 #####

#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps1014invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#101 taxa here in HI out of the total 125
psloinvaders <- subset_samples(ps1014invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#25 taxa here in LO out of the total 125. 
#Meaning they overlap in 1 ASV. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre110)))
hitax <- as.list(rownames(otu_table(pshipre110)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% lotax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#9 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% hitax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
# 11 taxa remaining 
#
ps1014sharedtax <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                  otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps1014sharedtax)$comparison <- "1014"
#
#      ##### D14 to D28 #####
#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps1428invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#7 taxa here in HI out of the total 32
psloinvaders <- subset_samples(ps1428invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#25 taxa here in LO out of the total 32. 
#Meaning they overlap in 0 ASVs. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre11014)))
hitax <- as.list(rownames(otu_table(pshipre11014)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% lotax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#2 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% hitax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
# 15 taxa remaining 
#
ps1428sharedtax <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                  otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps1428sharedtax)$comparison <- "1428"
#
#      ##### D28 to D41 #####
#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps2841invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#171 taxa here in HI out of the total 203
psloinvaders <- subset_samples(ps2841invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#39 taxa here in LO out of the total 203. 
#Meaning they overlap in 7 ASVs. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre1101428)))
hitax <- as.list(rownames(otu_table(pshipre1101428)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% lotax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#50 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% hitax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
# 16 taxa remaining 
#

ps2841sharedtax <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                  otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps2841sharedtax)$comparison <- "2841"
#
#
#         ##### Out of curiosity, how many of the 2841 invading HI ASVs were shared with previous HI timepoints? #####
# This time interval had a huge increase in worm-treated HI richness. 
# Remove all of the HI ASVs from previous timepoints from these invaders. 
tmpc <- subset(otu_table(pshiinvaders), !rownames(otu_table(pshiinvaders)) %in% hitax)
tmpc1 <- merge_phyloseq(tmpc, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#This is only 6 ASVs remaining out of the total 171. Meaning only 6 of the invaders had never been 
# in any previous HI samples. 


#      ##### Concatenate Invaders shared from the other microb into a large phyloseq object and df. #####
#into a phyloseq object
psallsharedtax <- merge_phyloseq(otu_table(pspre1sharedtax), tax_table(pspre1sharedtax), sample_data(pspre1sharedtax), 
                            otu_table(ps110sharedtax), tax_table(ps110sharedtax), sample_data(ps110sharedtax), 
                            otu_table(ps1014sharedtax), tax_table(ps1014sharedtax), sample_data(ps1014sharedtax), 
                            otu_table(ps1428sharedtax), tax_table(ps1428sharedtax), sample_data(ps1428sharedtax), 
                            otu_table(ps2841sharedtax), tax_table(ps2841sharedtax), sample_data(ps2841sharedtax)) 
sample_data(psallsharedtax)$readsShared <- sample_sums(psallsharedtax)
tmp <- estimate_richness(psallsharedtax, measures="Observed")
sample_data(psallsharedtax)$richnessShared <- tmp$Observed
sample_data(psallsharedtax)$mouse_comp <- paste(sample_data(psallsharedtax)$mouse, sample_data(psallsharedtax)$comparison, sep="_")

#Now pull it into a dataframe
dftmp <- data.frame(sample_data(psallsharedtax))
#Pull only the columns I want for plotting and future analysis. Going to merge this with dfinvaderdat. 
dfsharedtaxdat <- dftmp[,c("richnessShared", "readsShared", "mouse_comp", "sample", "day_num", "MT", "MTD")]

#


#   ##### Plot average richness of sharedtax over time by group #####
dfhitreat <- dfsharedtaxdat[dfsharedtaxdat$MT=="HI_worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(richnessShared)
dfhictrl <- dfsharedtaxdat[dfsharedtaxdat$MT=="HI_control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(richnessShared)
dflotreat <- dfsharedtaxdat[dfsharedtaxdat$MT=="LO_worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(richnessShared)
dfloctrl <- dfsharedtaxdat[dfsharedtaxdat$MT=="LO_control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(richnessShared)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessShared = mean(richnessShared, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessShared = mean(richnessShared, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessShared = mean(richnessShared, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessShared = mean(richnessShared, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessShared, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessShared, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessShared, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessShared, na.rm=TRUE))

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

dfhitreat.mean$SEH <- dfhitreat.mean$meanrichnessShared + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanrichnessShared + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanrichnessShared + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanrichnessShared + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanrichnessShared - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanrichnessShared - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanrichnessShared - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanrichnessShared - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as ASVrichnessShared-specific dfs. 
dfhitreat.mean.ASVrichnessShared <- dfhitreat.mean
##  day_num meanrichnessShared    SD     N   SEH    SEL
##1       1              0.444 0.616    18  0.590  0.299
##2      10              0.889 0.900    18  1.10   0.677
##3      14              0.467 0.640    15  0.632  0.301
##4      28              0.4   0.548     5  0.645  0.155
##5      41             13.8   6.14      8 15.9   11.6  
dfhictrl.mean.ASVrichnessShared <- dfhictrl.mean 
##1       1               3.5   2.38     4  4.69 2.31 
##2      10               6.25  1.5      4  7    5.5  #
##3      14               1     1.41     4  1.71 0.293
dflotreat.mean.ASVrichnessShared <- dflotreat.mean
##1       1              0.844 1.67     32  1.14   0.549
##2      10             17.2   7.78     32 18.5   15.8  
##3      14              0.312 0.535    32  0.407  0.218
##4      28              0.267 0.458    15  0.385  0.148
##5      41              1.52  1.12     23  1.76   1.29
dfloctrl.mean.ASVrichnessShared <- dfloctrl.mean 
##1       1              3      2.45     6  4    2    
##2      10              2.17   2.32     6  3.11 1.22 
##3      14              0.833  1.33     6  1.38 0.291
##4      28              3      3.35     6  4.37 1.63 
##5      41              1.67   1.86     6  2.43 0.907

#         ##### Plot as average over time with error bars #####
ASVrichnessShared <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.ASVrichnessShared$meanrichnessShared, x=dfhitreat.mean.ASVrichnessShared$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.ASVrichnessShared$meanrichnessShared, x=dfhictrl.mean.ASVrichnessShared$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.ASVrichnessShared$SEL, ymax=dfhitreat.mean.ASVrichnessShared$SEH, y=dfhitreat.mean.ASVrichnessShared$meanrichnessShared, x=dfhitreat.mean.ASVrichnessShared$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.ASVrichnessShared$SEL, ymax=dfhictrl.mean.ASVrichnessShared$SEH, y=dfhictrl.mean.ASVrichnessShared$meanrichnessShared, x=dfhictrl.mean.ASVrichnessShared$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.ASVrichnessShared$meanrichnessShared, x=dfhitreat.mean.ASVrichnessShared$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.ASVrichnessShared$meanrichnessShared, x=dfhictrl.mean.ASVrichnessShared$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.ASVrichnessShared$meanrichnessShared, x=dflotreat.mean.ASVrichnessShared$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.ASVrichnessShared$meanrichnessShared, x=dfloctrl.mean.ASVrichnessShared$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.ASVrichnessShared$SEL, ymax=dflotreat.mean.ASVrichnessShared$SEH, y=dflotreat.mean.ASVrichnessShared$meanrichnessShared, x=dflotreat.mean.ASVrichnessShared$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.ASVrichnessShared$SEL, ymax=dfloctrl.mean.ASVrichnessShared$SEH, y=dfloctrl.mean.ASVrichnessShared$meanrichnessShared, x=dfloctrl.mean.ASVrichnessShared$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.ASVrichnessShared$meanrichnessShared, x=dflotreat.mean.ASVrichnessShared$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.ASVrichnessShared$meanrichnessShared, x=dfloctrl.mean.ASVrichnessShared$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("ASV richness") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average Invader richness Shared w opp microb") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

ASVrichnessShared

#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_AverageASVrichnessShared_Resistance_lineplot_woD31.pdf", ASVrichnessShared, units="in", height=5, width=6)
#

#   ##### Plot average richness of sharedtax over time by group #####
dfhitreat <- dfsharedtaxdat[dfsharedtaxdat$MT=="HI_worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(richnessShared)
dfhictrl <- dfsharedtaxdat[dfsharedtaxdat$MT=="HI_control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(richnessShared)
dflotreat <- dfsharedtaxdat[dfsharedtaxdat$MT=="LO_worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(richnessShared)
dfloctrl <- dfsharedtaxdat[dfsharedtaxdat$MT=="LO_control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(richnessShared)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsShared = mean(readsShared, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsShared = mean(readsShared, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsShared = mean(readsShared, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanreadsShared = mean(readsShared, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsShared, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsShared, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsShared, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(readsShared, na.rm=TRUE))

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

dfhitreat.mean$SEH <- dfhitreat.mean$meanreadsShared + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanreadsShared + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanreadsShared + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanreadsShared + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanreadsShared - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanreadsShared - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanreadsShared - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanreadsShared - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as ASVreadsShared-specific dfs. 
dfhitreat.mean.ASVreadsShared <- dfhitreat.mean
## day_num meanrichness    SD     N   SEH    SEL
##1       1            2.72  6.25    18  4.19  1.25 
##2      10            1.44  1.79    18  1.87  1.02 
##3      14            1.13  1.92    15  1.63  0.637
##4      28            2.6   3.97     5  4.38  0.822
##5      41           45.6  27.0      8 55.2  36.1  
dfhictrl.mean.ASVreadsShared <- dfhictrl.mean 
##1       1             8.5  8.89     4 12.9  4.06 
##2      10            11.2  7.09     4 14.8  7.71 
##3      14             1.5  2.38     4  2.69 0.310
dflotreat.mean.ASVreadsShared <- dflotreat.mean
##1       1           1.81    4.97    32   2.69    0.934
##2      10           284.    202.      32 320.    248.   
##3      14           0.719   1.49    32   0.981   0.456
##4      28           2.33    5.84    15   3.84    0.826
##5      41           9.22    9.50    23  11.2     7.24 
dfloctrl.mean.ASVreadsShared <- dfloctrl.mean 
##1       1            7     8.65     6 10.5   3.47
##2      10           10.7  14.5      6 16.6   4.75
##3      14            3.33  5.20     6  5.46  1.21
##4      28            7.67 10.1      6 11.8   3.55
##5      41            3     3.74     6  4.53  1.47


#         ##### Plot as average over time with error bars #####
preadsShared <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.ASVreadsShared$meanreadsShared, x=dfhitreat.mean.ASVreadsShared$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.ASVreadsShared$meanreadsShared, x=dfhictrl.mean.ASVreadsShared$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.ASVreadsShared$SEL, ymax=dfhitreat.mean.ASVreadsShared$SEH, y=dfhitreat.mean.ASVreadsShared$meanreadsShared, x=dfhitreat.mean.ASVreadsShared$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.ASVreadsShared$SEL, ymax=dfhictrl.mean.ASVreadsShared$SEH, y=dfhictrl.mean.ASVreadsShared$meanreadsShared, x=dfhictrl.mean.ASVreadsShared$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.ASVreadsShared$meanreadsShared, x=dfhitreat.mean.ASVreadsShared$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.ASVreadsShared$meanreadsShared, x=dfhictrl.mean.ASVreadsShared$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.ASVreadsShared$meanreadsShared, x=dflotreat.mean.ASVreadsShared$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.ASVreadsShared$meanreadsShared, x=dfloctrl.mean.ASVreadsShared$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.ASVreadsShared$SEL, ymax=dflotreat.mean.ASVreadsShared$SEH, y=dflotreat.mean.ASVreadsShared$meanreadsShared, x=dflotreat.mean.ASVreadsShared$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.ASVreadsShared$SEL, ymax=dfloctrl.mean.ASVreadsShared$SEH, y=dfloctrl.mean.ASVreadsShared$meanreadsShared, x=dfloctrl.mean.ASVreadsShared$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.ASVreadsShared$meanreadsShared, x=dflotreat.mean.ASVreadsShared$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.ASVreadsShared$meanreadsShared, x=dfloctrl.mean.ASVreadsShared$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("Number readsShared (out of 3000)") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average readsShared") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

preadsShared

#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_AverageReadsShared_Resistance_lineplot_woD31.pdf", ASVreadsShared, units="in", height=5, width=6)
#

##### FINDING WHETHER INVADERS HAD BEEN IN THE SAME GROUP PREVIOUSLY #####
#      ##### D1 to D10 #####

#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps110invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#148 taxa here in HI out of the total 261
psloinvaders <- subset_samples(ps110invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#137 taxa here in LO out of the total 261. 
#Meaning they overlap in 24 ASVs. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(psprelo)))
hitax <- as.list(rownames(otu_table(psprehi)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% hitax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#108 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% lotax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
#16  taxa remaining 
#
ps110Reinvading <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                  otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps110Reinvading)$comparison <- "110"
#
#
#      ##### D10 to D14 #####

#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps1014invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#101 taxa here in HI out of the total 125
psloinvaders <- subset_samples(ps1014invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#25 taxa here in LO out of the total 125. 
#Meaning they overlap in 1 ASV. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre1)))
hitax <- as.list(rownames(otu_table(pshipre1)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% hitax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#88 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% lotax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
#15 taxa remaining 
#
ps1014Reinvading <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                   otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps1014Reinvading)$comparison <- "1014"
#
#      ##### D14 to D28 #####
#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps1428invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#7 taxa here in HI out of the total 32
psloinvaders <- subset_samples(ps1428invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#25 taxa here in LO out of the total 32. 
#Meaning they overlap in 0 ASVs. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre110)))
hitax <- as.list(rownames(otu_table(pshipre110)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% hitax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#4 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% lotax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
# 10 taxa remaining 
#
ps1428Reinvading <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                   otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps1428Reinvading)$comparison <- "1428"
#
#      ##### D28 to D41 #####
#Subset to hi and lo invaders for next portion of processing. 
pshiinvaders <- subset_samples(ps2841invaders, microb == "HI")
pshiinvaders <- prune_taxa(taxa_sums(pshiinvaders)>=1, pshiinvaders)
#171 taxa here in HI out of the total 203
psloinvaders <- subset_samples(ps2841invaders, microb == "LO")
psloinvaders <- prune_taxa(taxa_sums(psloinvaders)>=1, psloinvaders)
#39 taxa here in LO out of the total 203. 
#Meaning they overlap in 7 ASVs. 

#Grab the relevant taxa list to subtract from the invader df. 
lotax <- as.list(rownames(otu_table(pslopre11014)))
hitax <- as.list(rownames(otu_table(pshipre11014)))

#Keep the ASVs that ARE shared with the other df. 
#For HI animals. 
tmpa <- subset(otu_table(pshiinvaders), rownames(otu_table(pshiinvaders)) %in% hitax)
tmpa1 <- merge_phyloseq(tmpa, taxonomy=tax_table(pshiinvaders), metadata=sample_data(pshiinvaders))
#165 taxa remaining

#And for LO animals. 
tmpb <- subset(otu_table(psloinvaders), rownames(otu_table(psloinvaders)) %in% lotax)
tmpb1 <- merge_phyloseq(tmpb, taxonomy=tax_table(psloinvaders), metadata=sample_data(psloinvaders))
# 31 taxa remaining 
#

ps2841Reinvading <- merge_phyloseq(otu_table(tmpa1), tax_table(tmpa1),sample_data(tmpa1),
                                   otu_table(tmpb1), tax_table(tmpb1), sample_data(tmpb1))
sample_data(ps2841Reinvading)$comparison <- "2841"
#
#


#      ##### Concatenate Reinvading microbes into a large phyloseq object and df. #####
#into a phyloseq object
psallReinvading <- merge_phyloseq(otu_table(ps110Reinvading), tax_table(ps110Reinvading), sample_data(ps110Reinvading), 
                                 otu_table(ps1014Reinvading), tax_table(ps1014Reinvading), sample_data(ps1014Reinvading), 
                                 otu_table(ps1428Reinvading), tax_table(ps1428Reinvading), sample_data(ps1428Reinvading), 
                                 otu_table(ps2841Reinvading), tax_table(ps2841Reinvading), sample_data(ps2841Reinvading)) 
sample_data(psallReinvading)$readsReinvading <- sample_sums(psallReinvading)
tmp <- estimate_richness(psallReinvading, measures="Observed")
sample_data(psallReinvading)$richnessReinvading <- tmp$Observed
sample_data(psallReinvading)$mouse_comp <- paste(sample_data(psallReinvading)$mouse, sample_data(psallReinvading)$comparison, sep="_")

#Now pull it into a dataframe
dftmp <- data.frame(sample_data(psallReinvading))
#Pull only the columns I want for plotting and future analysis. Going to merge this with dfinvaderdat. 
dfReinvadingdat <- dftmp[,c("richnessReinvading", "readsReinvading", "mouse_comp", "sample", "day_num", "MT", "MTD")]

#



#   ##### Plot average richness of Reinvading over time by group #####
dfhitreat <- dfReinvadingdat[dfReinvadingdat$MT=="HI_worm",]
dfhitreat <- dfhitreat %>% tidyr::drop_na(richnessReinvading)
dfhictrl <- dfReinvadingdat[dfReinvadingdat$MT=="HI_control",]
dfhictrl <- dfhictrl %>% tidyr::drop_na(richnessReinvading)
dflotreat <- dfReinvadingdat[dfReinvadingdat$MT=="LO_worm",]
dflotreat <- dflotreat %>% tidyr::drop_na(richnessReinvading)
dfloctrl <- dfReinvadingdat[dfReinvadingdat$MT=="LO_control",]
dfloctrl <- dfloctrl %>% tidyr::drop_na(richnessReinvading)


#Calculate the means of each group. 
dfhitreat.mean <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessReinvading = mean(richnessReinvading, na.rm=TRUE))
dfhictrl.mean <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessReinvading = mean(richnessReinvading, na.rm=TRUE))
dflotreat.mean <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessReinvading = mean(richnessReinvading, na.rm=TRUE))
dfloctrl.mean <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(meanrichnessReinvading = mean(richnessReinvading, na.rm=TRUE))

#Get the standard deviation
dfhitreat.sd <- dfhitreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessReinvading, na.rm=TRUE))
dfhictrl.sd <- dfhictrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessReinvading, na.rm=TRUE))
dflotreat.sd <- dflotreat %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessReinvading, na.rm=TRUE))
dfloctrl.sd <- dfloctrl %>%
  group_by(day_num) %>%
  dplyr::summarize(SD = sd(richnessReinvading, na.rm=TRUE))

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

dfhitreat.mean$SEH <- dfhitreat.mean$meanrichnessReinvading + (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEH <- dfhictrl.mean$meanrichnessReinvading + (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEH <- dflotreat.mean$meanrichnessReinvading + (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEH <- dfloctrl.mean$meanrichnessReinvading + (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))

dfhitreat.mean$SEL <- dfhitreat.mean$meanrichnessReinvading - (dfhitreat.mean$SD/(sqrt(dfhitreat.mean$N)))
dfhictrl.mean$SEL <- dfhictrl.mean$meanrichnessReinvading - (dfhictrl.mean$SD/(sqrt(dfhictrl.mean$N)))
dflotreat.mean$SEL <- dflotreat.mean$meanrichnessReinvading - (dflotreat.mean$SD/(sqrt(dflotreat.mean$N)))
dfloctrl.mean$SEL <- dfloctrl.mean$meanrichnessReinvading - (dfloctrl.mean$SD/(sqrt(dfloctrl.mean$N)))


#save these as ASVrichnessReinvading-specific dfs. 
dfhitreat.mean.ASVrichnessReinvading <- dfhitreat.mean
##  day_num meanrichnessReinvading    SD     N   SEH    SEL
##1      10                   2.56  2.15    18  3.06  2.05 
##2      14                   5.6   3.64    15  6.54  4.66 
##3      28                   0.8   1.30     5  1.38  0.217
##4      41                  39.5  21.4      8 47.1  31.9  
dfhictrl.mean.ASVrichnessReinvading <- dfhictrl.mean 
##1      10                     35  11.5     4  40.8  29.2
##2      14                     13  10       4  18     8  
dflotreat.mean.ASVrichnessReinvading <- dflotreat.mean
##1      10                  0.688 0.931    32 0.852 0.523 
##2      14                  0.188 0.471    32 0.271 0.104 
##3      28                  0.2   0.414    15 0.307 0.0931
##4      41                  3.22  1.54     23 3.54  2.90  
dfloctrl.mean.ASVrichnessReinvading <- dfloctrl.mean 
##1      10                   2     1.41     6  2.58  1.42
##2      14                   2.33  2.58     6  3.39  1.28
##3      28                   2.17  1.72     6  2.87  1.46
##4      41                   4.17  2.04     6  5     3.33

#      ##### Plot as average over time with error bars #####
ASVrichnessReinvading <- ggplot() + 
  geom_vline(xintercept=0, color="darkgrey", linetype='solid') + 
  geom_vline(xintercept=24, color="darkgrey", linetype='solid') +
  theme_bw() + geom_line(aes(y=dfhitreat.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfhitreat.mean.ASVrichnessReinvading$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#FDD964") +
  geom_line(aes(y=dfhictrl.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfhictrl.mean.ASVrichnessReinvading$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#FDD964") + 
  geom_errorbar(aes(ymin=dfhitreat.mean.ASVrichnessReinvading$SEL, ymax=dfhitreat.mean.ASVrichnessReinvading$SEH, y=dfhitreat.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfhitreat.mean.ASVrichnessReinvading$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_errorbar(aes(ymin=dfhictrl.mean.ASVrichnessReinvading$SEL, ymax=dfhictrl.mean.ASVrichnessReinvading$SEH, y=dfhictrl.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfhictrl.mean.ASVrichnessReinvading$day_num), alpha=1, color="#FDD964", width=0.8) + 
  geom_point(aes(y=dfhitreat.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfhitreat.mean.ASVrichnessReinvading$day_num),shape=15, alpha=1, size=3, na.rm = TRUE, color="#FDD964") + 
  geom_point(aes(y=dfhictrl.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfhictrl.mean.ASVrichnessReinvading$day_num), shape=22, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#FDD964") + 
  geom_line(aes(y=dflotreat.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dflotreat.mean.ASVrichnessReinvading$day_num),alpha=0.9, linewidth=0.75, na.rm = TRUE, color="#C66EF5") +
  geom_line(aes(y=dfloctrl.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfloctrl.mean.ASVrichnessReinvading$day_num), alpha=0.9, linewidth=0.6, linetype=2,na.rm = TRUE, color="#C66EF5") + 
  geom_errorbar(aes(ymin=dflotreat.mean.ASVrichnessReinvading$SEL, ymax=dflotreat.mean.ASVrichnessReinvading$SEH, y=dflotreat.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dflotreat.mean.ASVrichnessReinvading$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_errorbar(aes(ymin=dfloctrl.mean.ASVrichnessReinvading$SEL, ymax=dfloctrl.mean.ASVrichnessReinvading$SEH, y=dfloctrl.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfloctrl.mean.ASVrichnessReinvading$day_num), alpha=1, color="#C66EF5", width=0.8) + 
  geom_point(aes(y=dflotreat.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dflotreat.mean.ASVrichnessReinvading$day_num),shape=16, alpha=1, size=3, na.rm = TRUE, color="#C66EF5") + 
  geom_point(aes(y=dfloctrl.mean.ASVrichnessReinvading$meanrichnessReinvading, x=dfloctrl.mean.ASVrichnessReinvading$day_num), shape=21, fill= '#ffffff',alpha=1, size=3,na.rm = TRUE, color="#C66EF5") + 
  ylab("ASV richness") +
  xlab("Sampling timepoint (dpi)") +
  ggtitle("Average richness ASVs Reinvading from past timepoints") +
  theme(plot.title=element_text(size=14, face="bold"), 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13), 
        legend.position="none")
#scale_x_continuous(name="Time in day_nums", limits=c(-22, 17), breaks=c(-21, -7, 0, 4, 8, 12, 16)) +
#ylim(150, 350) + 
#annotate(geom="text", x=9, y=155, label="Experiment", color="purple", size=5) + 
#annotate(geom="text", x=-11.5, y=155, label="Acclimation", color="purple", size=5) 

ASVrichnessReinvading

#ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/psexp3000_AverageASVrichnessReinvading_Resistance_lineplot_woD31.pdf", ASVrichnessReinvading, units="in", height=5, width=6)
#



##### CREATING A MASTER DF FOR ALL INVADING, LOST, SHARED ASVs #####
#Finally, pull the reads and richness data from each of the "taxa lost" and "shared with 
# other microb type" and also create another column for "not shared with the other microbiome type".
# Need to pull out the reads and richness from each of the phyloseq objects, make unique names
# for the columns from each type of analysis. 

#column names are based on reads and richness, with "sample" being the merge column. Also want to pull treatment, microb, MT, day, day_num, mouse from the first one. 

#First, create df's and bind together using the invader dfs as the starting point. 
#Already made the invader dataframe up in the first section. Just need to add columns to it. 

colnames(dfinvaderdat)
colnames(dflostdat)
colnames(dfsharedtaxdat)
colnames(dfReinvadingdat)
#Combine the df's in 2 steps. 
dftmp <- merge(dfinvaderdat, dfsharedtaxdat[,c(1:4)], by="sample", all.x=TRUE)
dftmp1 <- merge(dftmp, dfReinvadingdat[,c(1,2,4)], by="sample", all.x=TRUE)
dfInvaderLostShared1 <- merge(dftmp1, dflostdat, by="mouse_comp", all=TRUE) #This has lots of NA's for because there are 19 samples that have values for lost dat but not the other 2 dfs. 
#Now create one that doesn't keep all the extra rows from dflostdat. 
dfInvaderLostShared2 <- merge(dftmp1, dflostdat, by="mouse_comp", all.x=TRUE) #This is only for samples that have invaderdat, dropping the lostdat if no other values. 
#

#Aggregate the data to easily get means and sd's by group. 
avg <- aggregate(richnessInvaders ~ MTD, data=dfInvaderLostShared2, FUN="mean")
sd <- aggregate(richnessInvaders ~ MTD, data=dfInvaderLostShared2, FUN="sd")
summaryInvaders <- cbind(avg, sd$richnessInvaders)
#View(summaryInvaders)

#For ASVs that may be sourced from the other microbiome type. 
avg <- aggregate(richnessShared ~ MTD, data=dfInvaderLostShared2, FUN="mean")
sd <- aggregate(richnessShared ~ MTD, data=dfInvaderLostShared2, FUN="sd")
summaryShared <- cbind(avg, sd$richnessShared)
#View(summaryShared)

#For lost ASVs
avg <- aggregate(richnessLost ~ MTD, data=dfInvaderLostShared1, FUN="mean")
sd <- aggregate(richnessLost ~ MTD, data=dfInvaderLostShared1, FUN="sd")
summaryLost <- cbind(avg, sd$richnessLost)
View(summaryLost)

#For Reinvading ASVs. 
avg <- aggregate(richnessReinvading ~ MTD, data=dfInvaderLostShared2, FUN="mean")
sd <- aggregate(richnessReinvading ~ MTD, data=dfInvaderLostShared2, FUN="sd")
summaryReinvading <- cbind(avg, sd$richnessReinvading)
#View(summaryReinvading)

#Now do this for raw richness so that it's easy to access
dfmetadatwhole$MTD <- paste(dfmetadatwhole$MT, dfmetadatwhole$day, sep="_")
avg <- aggregate(Observed ~ MTD, data=dfmetadatwhole, FUN="mean")
sd <- aggregate(Observed ~ MTD, data=dfmetadatwhole, FUN="sd")
summaryRawRichness <- cbind(avg, sd$Observed)
View(summaryRawRichness)

#How about the number of reads that can be attributed to the original 11 PhyLo11B taxa?
dfrare11 <- data.frame(sample_data(rare11), 
                     estimate_richness(rare11, measures="Observed"), 
                     "reads"=sample_sums(rare11))
avg <- aggregate(reads ~ MTD, data=dfrare11, FUN="mean")
sd <- aggregate(reads ~ MTD, data=dfrare11, FUN="sd")
summaryPhyLo11BReads <- cbind(avg, sd$reads)
View(summaryPhyLo11BReads)


View(avg)

#   ##### Plotting from the summary data #####

#      ##### Use the summarized data to plot barplots of average richness of invaders for each group #####
#I updated the summarized Invaders spreadsheet in Excel to add metadata columns. 
#Now reading it back in. 
summaryInvaders <- read.csv("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/summaryInvadersUpdated.csv")
#
mycols <- c("LO_worm"= "#C66EF5",
            "LO_control" = "#C66EF5", 
            "HI_worm" = "#FDD964",
            "HI_control" = "#FDD964" )

myfill <- c('LO_worm' = "#C66EF5", 
            'LO_control' = "#FFFFFF", 
            'HI_worm' = "#FDD964",
            'HI_control' = "#FFFFFF")
#

summaryInvaders$MTD <- factor(summaryInvaders$MTD, 
                              levels=c("LO_worm_D1",  "LO_worm_D10" ,   "LO_worm_D14" , "LO_worm_D28" , "LO_worm_D41",
                                       "LO_control_D1", "LO_control_D10" , "LO_control_D14","LO_control_D28" ,"LO_control_D41",    
                                       "HI_worm_D1" , "HI_worm_D10", "HI_worm_D14",   "HI_worm_D28",   "HI_worm_D41", 
                                       "HI_control_D1", "HI_control_D10" ,"HI_control_D14"))

pinvaderbarpre1 <- ggplot(summaryInvaders[summaryInvaders$day=="D1",], aes(x=MTD, y=richnessInvaders)) + 
  geom_bar(stat="identity", aes(color=MT, fill=MT)) +
  geom_errorbar(aes(ymin = richnessInvaders-sd, ymax= richnessInvaders+sd, color=MT), width=0.5) + 
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=myfill)+
  ylim(0,55) +
  ylab("Number of invading ASVs") +
  xlab("Group") +
  ggtitle("Number of invaders since previous sample") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
# facet_wrap(~day, scales="free_x", ncol=5) 
pinvaderbarpre1

pinvaderbar110 <- ggplot(summaryInvaders[summaryInvaders$day=="D10",], aes(x=MTD, y=richnessInvaders)) + 
  geom_bar(stat="identity", aes(color=MT, fill=MT)) +
  geom_errorbar(aes(ymin = richnessInvaders-sd, ymax= richnessInvaders+sd, color=MT), width=0.5) + 
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=myfill)+
  ylim(0,55) +
  ylab("Number of invading ASVs") +
  xlab("Group") +
  ggtitle("Number of invaders since previous sample") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
pinvaderbar110

pinvaderbar1014 <- ggplot(summaryInvaders[summaryInvaders$day=="D14",], aes(x=MTD, y=richnessInvaders)) + 
  geom_bar(stat="identity", aes(color=MT, fill=MT)) +
  geom_errorbar(aes(ymin = richnessInvaders-sd, ymax= richnessInvaders+sd, color=MT), width=0.5) + 
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=myfill)+
  ylim(0,55) +
  ylab("Number of invading ASVs") +
  xlab("Group") +
  ggtitle("Number of invaders since previous sample") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
pinvaderbar1014

pinvaderbar1428 <- ggplot(summaryInvaders[summaryInvaders$day=="D28",], aes(x=MTD, y=richnessInvaders)) + 
  geom_bar(stat="identity", aes(color=MT, fill=MT)) +
  geom_errorbar(aes(ymin = richnessInvaders-sd, ymax= richnessInvaders+sd, color=MT), width=0.5) + 
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=myfill)+
  ylim(0,10) +
  ylab("Number of invading ASVs") +
  xlab("Group") +
  ggtitle("Number of invaders since previous sample") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
pinvaderbar1428

pinvaderbar2841 <- ggplot(summaryInvaders[summaryInvaders$day=="D41",], aes(x=MTD, y=richnessInvaders)) + 
  geom_bar(stat="identity", aes(color=MT, fill=MT)) +
  geom_errorbar(aes(ymin = richnessInvaders-sd, ymax= richnessInvaders+sd, color=MT), width=0.5) + 
  scale_colour_manual(values=mycols) +
  scale_fill_manual(values=myfill)+
  ylim(0,65) +
  ylab("Number of invading ASVs") +
  xlab("Group") +
  ggtitle("Number of invaders since previous sample") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
pinvaderbar2841


#Save the faceted plot
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/InvaderBarplotpre1.pdf", pinvaderbarpre1, units="in", height=3, width=2)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/InvaderBarplot110.pdf", pinvaderbar110, units="in", height=3, width=2)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/InvaderBarplot1014.pdf", pinvaderbar1014, units="in", height=3, width=2)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/InvaderBarplot1428.pdf", pinvaderbar1428, units="in", height=3, width=1.5)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resistance/InvaderBarplot2841.pdf", pinvaderbar2841, units="in", height=3, width=1.5)
#

#
#      ##### Now creating a dot plot of the numbers of diff abund with invasion info included #####
View(ResistDiffAbundBC)

mycols <- c("LO_worm"= "#C66EF5",
            "LO_sham" = "#C66EF5", 
            "HI_worm" = "#FDD964",
            "HI_sham" = "#FDD964" )
#
myfill <- c('LO_worm_no' = "#C66EF5", 
             'LO_sham_no' = "#FFFFFF", 
             'HI_worm_no' = "#FDD964",
             'HI_sham_no' = "#FFFFFF", 
             'LO_worm_yes' = "#000000", 
             'LO_sham_yes' = "#000000", 
             'HI_worm_yes' = "#000000",
             'HI_sham_yes' = "#000000")
#
#Write out the differential abundance df and update so that there are blank rows for times that don't have anything in them. 
#write.csv(ResistDiffAbundBC, "~/Documents/1U4U_16S/GG2_reclassification/DESeq2/ResistDiffAbundBCPlusInvaders.csv")
#Add the extra rows as placeholders for MTD's that don't have diff abund ASVs. 
ResistDiffAbundBCtoplot <- read.csv("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/ResistDiffAbundBCPlusInvaders.csv")

# Create a new column for labeling purposes. 
ResistDiffAbundBCtoplot$label <- paste(ResistDiffAbundBCtoplot$MT, ResistDiffAbundBCtoplot$invading, sep = "_")

ResistDiffAbundBCtoplot$MTD <- factor(ResistDiffAbundBCtoplot$MTD, 
                                levels=c('LO_worm_pre1', 'LO_worm_110', 'LO_worm_1014', 'LO_worm_1428', 'LO_worm_2841', 
                                         'LO_sham_pre1', 'LO_sham_110', 'LO_sham_1014', 'LO_sham_1428', 'LO_sham_2841', 
                                         'HI_worm_pre1', 'HI_worm_110' , 'HI_worm_1014' , 'HI_worm_1428', 'HI_worm_2841', 
                                         'HI_sham_pre1', 'HI_sham_110', 'HI_sham_1014'))
#
pinvaderpointpre1 <- ggplot(ResistDiffAbundBCtoplot[ResistDiffAbundBCtoplot$comparison=="pre1",], aes(x=MTD, y=log2FoldChange)) +
  geom_hline(yintercept=0, color="black", linetype='dashed') + 
  geom_beeswarm(aes(color=MT, fill=label, shape= MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(22, 22, 21, 21)) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values=myfill) +
  ylim(-10,25) +
  ylab("log2 fold change") +
  xlab("Group") +
  ggtitle("Differential abundance with invaders") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
  #facet_wrap(~comparison, scales="free_x", ncol=5)
pinvaderpointpre1

#
pinvaderpoint110 <- ggplot(ResistDiffAbundBCtoplot[ResistDiffAbundBCtoplot$comparison=="110",], aes(x=MTD, y=log2FoldChange)) +
  geom_hline(yintercept=0, color="black", linetype='dashed') + 
  geom_beeswarm(aes(color=MT, fill=label, shape= MT), cex=2.5, size=3, alpha=0.9) +
  scale_shape_manual(values=c(22, 22, 21, 21)) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values=myfill) +
  ylim(-10,25) +
  ylab("log2 fold change") +
  xlab("Group") +
  ggtitle("Differential abundance with invaders") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
#facet_wrap(~comparison, scales="free_x", ncol=5)
pinvaderpoint110

#
pinvaderpoint1014 <- ggplot(ResistDiffAbundBCtoplot[ResistDiffAbundBCtoplot$comparison=="1014",], aes(x=MTD, y=log2FoldChange)) +
  geom_hline(yintercept=0, color="black", linetype='dashed') + 
  geom_beeswarm(aes(color=MT, fill=label, shape= MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(22, 22, 21, 21)) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values=myfill) +
  ylim(-10,25) +
  ylab("log2 fold change") +
  xlab("Group") +
  ggtitle("Differential abundance with invaders") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
#facet_wrap(~comparison, scales="free_x", ncol=5)
pinvaderpoint1014

pinvaderpoint1428 <- ggplot(ResistDiffAbundBCtoplot[ResistDiffAbundBCtoplot$comparison=="1428",], aes(x=MTD, y=log2FoldChange)) +
  geom_hline(yintercept=0, color="black", linetype='dashed') + 
  geom_beeswarm(aes(color=MT, fill=label, shape= MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(22, 21, 21)) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values=myfill) +
  ylim(-10,25) +
  ylab("log2 fold change") +
  xlab("Group") +
  ggtitle("Differential abundance with invaders") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
#facet_wrap(~comparison, scales="free_x", ncol=5)
pinvaderpoint1428

pinvaderpoint2841 <- ggplot(ResistDiffAbundBCtoplot[ResistDiffAbundBCtoplot$comparison=="2841",], aes(x=MTD, y=log2FoldChange)) +
  geom_hline(yintercept=0, color="black", linetype='dashed') + 
  geom_beeswarm(aes(color=MT, fill=label, shape= MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(22, 21, 21)) +
  scale_color_manual(values=mycols) +
  scale_fill_manual(values=myfill) +
  ylim(-10,25) +
  ylab("log2 fold change") +
  xlab("Group") +
  ggtitle("Differential abundance with invaders") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  ggpubr::rotate_x_text(40) 
#facet_wrap(~comparison, scales="free_x", ncol=5)
pinvaderpoint2841

#
#Save the faceted plot
ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/InvaderDiffAbundPointpre1.pdf", pinvaderpointpre1, units="in", height=6, width=4)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/InvaderDiffAbundPoint110.pdf", pinvaderpoint110, units="in", height=6, width=4)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/InvaderDiffAbundPoint1014.pdf", pinvaderpoint1014, units="in", height=6, width=4)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/InvaderDiffAbundPoint1428.pdf", pinvaderpoint1428, units="in", height=6, width=3)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/DESeq2/InvaderDiffAbundPoint2841.pdf", pinvaderpoint2841, units="in", height=6, width=3)
#

#
#   ##### Stats on the master df #####
#      ##### Running wilcoxons for reads of invading taxa.  #####
#Subset the main df for each day. 
dfD1invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D1",]
dfD10invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D10",]
dfD14invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D14",]
dfD28invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D28",]
dfD41invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D41",]


#Testing change in ASV counts between pretreatment and D1. 
df1 <- dfD1invaders[dfD1invaders$treatment == "worm",]
wilcox.test(readsInvaders ~ microb, data=df1)
#W = 482, p-value = 5.767e-05
df2 <- dfD1invaders[dfD1invaders$microb == "LO",]
wilcox.test(readsInvaders ~ treatment, data=df2)
#W = 152.5, p-value = 0.01807
df3 <- dfD1invaders[dfD1invaders$microb == "HI",]
wilcox.test(readsInvaders ~ treatment, data=df3)
#W = 62, p-value = 0.02917
ps <- c(0.00005767 , 0.01807, 0.02917)
p.adjust(ps, "holm")
## 0.00017301 0.03614000 0.03614000

#Testing change in ASV counts between D1 and D10. 
df1 <- dfD10invaders[dfD10invaders$treatment == "worm",]
wilcox.test(readsInvaders ~ microb, data=df1)
#W = 0, p-value = 6.168e-09
df2 <- dfD10invaders[dfD10invaders$microb == "LO",]
wilcox.test(readsInvaders ~ treatment, data=df2)
#W = 0, p-value = 7.245e-07
df3 <- dfD10invaders[dfD10invaders$microb == "HI",]
wilcox.test(readsInvaders ~ treatment, data=df3)
#W = 72, p-value = 0.002441
ps <- c(0.000000006168, 0.0000007245, 0.002441)
p.adjust(ps, "holm")
## 1.8504e-08 1.4490e-06 2.4410e-03


#Testing change in ASV counts between D10 and D14. 
df1 <- dfD14invaders[dfD14invaders$treatment == "worm",]
wilcox.test(readsInvaders ~ microb, data=df1)
#W = 458.5, p-value = 1.487e-07
df2 <- dfD14invaders[dfD14invaders$microb == "LO",]
wilcox.test(readsInvaders ~ treatment, data=df2)
#W = 151, p-value = 0.01325
df3 <- dfD14invaders[dfD14invaders$microb == "HI",]
wilcox.test(readsInvaders ~ treatment, data=df3)
#W = 55.5, p-value = 0.01223
ps <- c(0.0000001487 , 0.01325, 0.01223)
p.adjust(ps, "holm")
## 4.461e-07 2.446e-02 2.446e-02

#Testing change in ASV counts between D14 and D28. 
df1 <- dfD28invaders[dfD28invaders$treatment == "worm",]
wilcox.test(readsInvaders ~ microb, data=df1)
#W = 39.5, p-value = 0.8934
df2 <- dfD28invaders[dfD28invaders$microb == "LO",]
wilcox.test(readsInvaders ~ treatment, data=df2)
#W = 52.5, p-value = 0.5801
#
ps <- c(0.8934 , 0.5801)
p.adjust(ps, "holm")
## 1 1

#Testing change in evenness between D28 and D41
df1 <- dfD41invaders[dfD41invaders$treatment == "worm",]
wilcox.test(readsInvaders ~ microb, data=df1)
#W = 180, p-value = 7.655e-05
df2 <- dfD41invaders[dfD41invaders$microb == "LO",]
wilcox.test(readsInvaders ~ treatment, data=df2)
#W = 63, p-value = 0.7667
ps <- c(0.00007655, 0.7667)
p.adjust(ps, "holm")
#0.0001531 0.7667000
#



#      ##### Running wilcoxons for richness of invading taxa.  #####
#Subset the main df for each day. 
dfD1invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D1",]
dfD10invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D10",]
dfD14invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D14",]
dfD28invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D28",]
dfD41invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D41",]


#Testing change in ASV counts between pretreatment and D1. 
df1 <- dfD1invaders[dfD1invaders$treatment == "worm",]
wilcox.test(richnessInvaders ~ microb, data=df1)
#W = 491.5, p-value = 2.261e-05
df2 <- dfD1invaders[dfD1invaders$microb == "LO",]
wilcox.test(richnessInvaders ~ treatment, data=df2)
#W = 160.5, p-value = 0.006481
df3 <- dfD1invaders[dfD1invaders$microb == "HI",]
wilcox.test(richnessInvaders ~ treatment, data=df3)
#W = 67.5, p-value = 0.007992
ps <- c(0.00002261 , 0.006481, 0.007992)
p.adjust(ps, "holm")
## 0.00006783 0.01296200 0.01296200

#Testing change in ASV counts between D1 and D10. 
df1 <- dfD10invaders[dfD10invaders$treatment == "worm",]
wilcox.test(richnessInvaders ~ microb, data=df1)
#W = 0, p-value = 5.985e-09
df2 <- dfD10invaders[dfD10invaders$microb == "LO",]
wilcox.test(richnessInvaders ~ treatment, data=df2)
#W = 0, p-value = 0.0001295
df3 <- dfD10invaders[dfD10invaders$microb == "HI",]
wilcox.test(richnessInvaders ~ treatment, data=df3)
#W = 72, p-value = 0.002332
ps <- c(0.000000005985, 0.0001295, 0.002332)
p.adjust(ps, "holm")
## 1.7955e-08 2.5900e-04 2.3320e-03


#Testing change in ASV counts between D10 and D14. 
df1 <- dfD14invaders[dfD14invaders$treatment == "worm",]
wilcox.test(richnessInvaders ~ microb, data=df1)
#W = 467, p-value = 4.584e-08
df2 <- dfD14invaders[dfD14invaders$microb == "LO",]
wilcox.test(richnessInvaders ~ treatment, data=df2)
#W = 152.5, p-value = 0.01069
df3 <- dfD14invaders[dfD14invaders$microb == "HI",]
wilcox.test(richnessInvaders ~ treatment, data=df3)
#W = 49, p-value = 0.06315
ps <- c(0.00000004584 , 0.01069, 0.06315)
p.adjust(ps, "holm")
## 1.3752e-07 2.1380e-02 6.3150e-02

#Testing change in ASV counts between D14 and D28. 
df1 <- dfD28invaders[dfD28invaders$treatment == "worm",]
wilcox.test(richnessInvaders ~ microb, data=df1)
#W = 38, p-value = 1
df2 <- dfD28invaders[dfD28invaders$microb == "LO",]
wilcox.test(richnessInvaders ~ treatment, data=df2)
#W = 69.5, p-value = 0.04866
#
ps <- c(1 , 0.04866)
p.adjust(ps, "holm")
## 1.00000 0.09732

#Testing change in evenness between D28 and D41
df1 <- dfD41invaders[dfD41invaders$treatment == "worm",]
wilcox.test(richnessInvaders ~ microb, data=df1)
#W = 184, p-value = 2.95e-05
df2 <- dfD41invaders[dfD41invaders$microb == "LO",]
wilcox.test(richnessInvaders ~ treatment, data=df2)
#W = 79.5, p-value = 0.5828
ps <- c(0.0000295 , 0.5828)
p.adjust(ps, "holm")
#0.000059 0.582800
#
#      ##### Running wilcoxons for reads of Taxa from opposite microbiome type.  #####
#Subset the main df for each day. 
dfD1invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D1",]
dfD10invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D10",]
dfD14invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D14",]
dfD28invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D28",]
dfD41invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D41",]


#Testing change in ASV counts between pretreatment and D1. 
df1 <- dfD1invaders[dfD1invaders$treatment == "worm",]
wilcox.test(readsShared ~ microb, data=df1)
#W = 294.5, p-value = 0.8907
df2 <- dfD1invaders[dfD1invaders$microb == "LO",]
wilcox.test(readsShared ~ treatment, data=df2)
#W = 156.5, p-value = 0.009635
df3 <- dfD1invaders[dfD1invaders$microb == "HI",]
wilcox.test(readsShared ~ treatment, data=df3)
#W = 60, p-value = 0.03204
ps <- c(0.8907 , 0.009635, 0.03204)
p.adjust(ps, "holm")
## 0.890700 0.028905 0.064080

#Testing change in ASV counts between D1 and D10. 
df1 <- dfD10invaders[dfD10invaders$treatment == "worm",]
wilcox.test(readsShared ~ microb, data=df1)
#W = 0, p-value = 5.822e-09
df2 <- dfD10invaders[dfD10invaders$microb == "LO",]
wilcox.test(readsShared ~ treatment, data=df2)
#W = 0, p-value = 0.0001316
df3 <- dfD10invaders[dfD10invaders$microb == "HI",]
wilcox.test(readsShared ~ treatment, data=df3)
#W = 70, p-value = 0.003518
ps <- c(0.000000005822, 0.0001316, 0.003518)
p.adjust(ps, "holm")
## 1.7466e-08 2.6320e-04 3.5180e-03


#Testing change in ASV counts between D10 and D14. 
df1 <- dfD14invaders[dfD14invaders$treatment == "worm",]
wilcox.test(readsShared ~ microb, data=df1)
#W = 270, p-value = 0.4145
df2 <- dfD14invaders[dfD14invaders$microb == "LO",]
wilcox.test(readsShared ~ treatment, data=df2)
#W = 110, p-value = 0.4994
df3 <- dfD14invaders[dfD14invaders$microb == "HI",]
wilcox.test(readsShared ~ treatment, data=df3)
#W = 32.5, p-value = 0.8233
ps <- c(0.4145 , 0.4994, 0.8233)
p.adjust(ps, "holm")
## 1.0000 1.0000 1.0000

#Testing change in ASV counts between D14 and D28. 
df1 <- dfD28invaders[dfD28invaders$treatment == "worm",]
wilcox.test(readsShared ~ microb, data=df1)
#W = 42.5, p-value = 0.628
df2 <- dfD28invaders[dfD28invaders$microb == "LO",]
wilcox.test(readsShared ~ treatment, data=df2)
#W = 65, p-value = 0.08222
#
ps <- c(0.628 , 0.08222)
p.adjust(ps, "holm")
## 0.62800 0.16444

#Testing change in evenness between D28 and D41
df1 <- dfD41invaders[dfD41invaders$treatment == "worm",]
wilcox.test(readsShared ~ microb, data=df1)
#W = 176, p-value = 0.0001601
df2 <- dfD41invaders[dfD41invaders$microb == "LO",]
wilcox.test(readsShared ~ treatment, data=df2)
#W = 36.5, p-value = 0.08346
ps <- c(0.0001601 , 0.08346)
p.adjust(ps, "holm")
#0.0003202 0.0834600
#


#      ##### Running wilcoxons for richness od taxa from opposite microbiome type.  #####
#Subset the main df for each day. 
dfD1invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D1",]
dfD10invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D10",]
dfD14invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D14",]
dfD28invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D28",]
dfD41invaders <- dfInvaderLostShared2[dfInvaderLostShared2$day=="D41",]


#Testing change in ASV counts between pretreatment and D1. 
df1 <- dfD1invaders[dfD1invaders$treatment == "worm",]
wilcox.test(richnessShared ~ microb, data=df1)
#W = 271, p-value = 0.7024
df2 <- dfD1invaders[dfD1invaders$microb == "LO",]
wilcox.test(richnessShared ~ treatment, data=df2)
#W = 161, p-value = 0.005117
df3 <- dfD1invaders[dfD1invaders$microb == "HI",]
wilcox.test(richnessShared ~ treatment, data=df3)
#W = 67.5, p-value = 0.004064
ps <- c(0.7024 , 0.005117, 0.004064)
p.adjust(ps, "holm")
## 0.702400 0.012192 0.012192

#Testing change in ASV counts between D1 and D10. 
df1 <- dfD10invaders[dfD10invaders$treatment == "worm",]
wilcox.test(richnessShared ~ microb, data=df1)
#W = 0, p-value = 5.447e-09
df2 <- dfD10invaders[dfD10invaders$microb == "LO",]
wilcox.test(richnessShared ~ treatment, data=df2)
#W = 0, p-value = 0.000127
df3 <- dfD10invaders[dfD10invaders$microb == "HI",]
wilcox.test(richnessShared ~ treatment, data=df3)
#W = 72, p-value = 0.001766
ps <- c(0.000000005447, 0.000127, 0.001766)
p.adjust(ps, "holm")
## 1.6341e-08 2.5400e-04 1.7660e-03


#Testing change in ASV counts between D10 and D14. 
df1 <- dfD14invaders[dfD14invaders$treatment == "worm",]
wilcox.test(richnessShared ~ microb, data=df1)
#W = 270, p-value = 0.4085
df2 <- dfD14invaders[dfD14invaders$microb == "LO",]
wilcox.test(richnessShared ~ treatment, data=df2)
#W = 109.5, p-value = 0.5128
df3 <- dfD14invaders[dfD14invaders$microb == "HI",]
wilcox.test(richnessShared ~ treatment, data=df3)
#W = 35.5, p-value = 0.5704
ps <- c(0.4085 , 0.5128, 0.5704)
p.adjust(ps, "holm")
## 1.0000 1.0000 1.0000

#Testing change in ASV counts between D14 and D28. 
df1 <- dfD28invaders[dfD28invaders$treatment == "worm",]
wilcox.test(richnessShared ~ microb, data=df1)
#W = 42.5, p-value = 0.6211
df2 <- dfD28invaders[dfD28invaders$microb == "LO",]
wilcox.test(richnessShared ~ treatment, data=df2)
#W = 71, p-value = 0.02252
#
ps <- c(0.6211 , 0.02252)
p.adjust(ps, "holm")
## 0.62110 0.04504

#Testing change in evenness between D28 and D41
df1 <- dfD41invaders[dfD41invaders$treatment == "worm",]
wilcox.test(richnessShared ~ microb, data=df1)
#W = 184, p-value = 2.356e-05
df2 <- dfD41invaders[dfD41invaders$microb == "LO",]
wilcox.test(richnessShared ~ treatment, data=df2)
#W = 63, p-value = 0.7563
ps <- c(0.00002356 , 0.7563)
p.adjust(ps, "holm")
#4.712e-05 7.563e-01
#




