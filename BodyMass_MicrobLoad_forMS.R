#body mass and microbial load analyses and visualizations

########### Body mass analysis and plotting #####
library(ggbeeswarm)
library(ggplot2)
library(nlme)
library(emmeans)


dat <- read.csv("mouse_body_mass.csv")
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
  labs(x="Sampling timepoint", y="Body mass (g)", title="Body Mass All Animals 1U4U")


ggsave("~/Documents/1U4U_16S/GG2_reclassification/bodymassBoxplots_psimptimepoints.pdf", units="in", height=5, width=9)


#modelling change in body mass over time in worm-treated animals, separately for LO and HI

mem.bod <- lme(mass ~ microb + treatment + day_num + microb*day_num, random=~1|mouse, data=dat)
mem.bod
anova(mem.bod)



datlo <- dat[dat$microb=="lo",]
dathi <- dat[dat$microb=="hi",]

mean(dathi$mass[dathi$day_num=="-5"],)
sd(dathi$mass[dathi$day_num=="-5"],)



mem.bod <- lme(mass ~ treatment + day, random=~1|mouse, data=datlo)
mem.bod
anova(mem.bod)


mem.bod <- lme(mass ~ treatment + day_num, random=~1|mouse, data=dathi)
mem.bod
anova(mem.bod)

#post-hoc testing with emmeans
emmeans(mem.bod, pairwise ~ day)



#


########### Microbial load analyses #####
#qPCR analysis, just to see whether there are any differences in microbial load. 
# Comparing the outcomes of microbial load by group. 

#He standardized the gDNA concentrations for all samples, and then ran qPCR on them. Apparently
# the qPCR machine has software that does some QC on its own, then can get the total 16S copy
# number in each sample. This is the column that we're most interested in. 

#Create your working directory
setwd("~/Documents/1U4U_qPCR")

#Read in the whole df. 

allqpcr <- read.csv("./1U4U_Full16S_metadata_ByOrigID_wBactQuant_wAll18684.csv") #original data from Zac. 
newqpcr <- read.csv("map_1U4U_Full16S_metadata_ByOrigID_wBactQuant.csv") 

#What are the variables that we care about?
colnames(allqpcr)
colnames(newqpcr)

#
# ##### Add a couple important columns #####

#Add read count to see whether overall sample quality is affecting qPCR outcomes. 
#Pull the read counts. 
dfps1 <- read.csv("ps1_metadata.csv")
div_day <- paste(dfps1$microb, sep="_", dfps1$day)
reads <- cbind("GnomexID"=dfps1$GNomexID, "ASVs"=dfps1$Observed,
               "count"=dfps1$sum, div_day)


#Merge the read counts with the allqpcr df. 
allqpcr.1 <- merge(allqpcr, reads, by="GnomexID")
allqpcr.1$count <- as.integer(allqpcr.1$count) #Make sure it's numeric
allqpcr.1$ASVs <- as.integer(allqpcr.1$ASVs) #Make sure it's numeric
allqpcr.1$microb <- factor(allqpcr.1$microb, levels=c("LO", "HI"))
head(allqpcr.1)
#Create a new column that is 16S copy number divided by number of reads to standardize? 
allqpcr.1$copiesbyreads <- allqpcr.1$copies_16S/allqpcr.1$count
max(allqpcr.1$copiesbyreads)
boxplot(allqpcr.1$copiesbyreads ~ allqpcr.1$day)



#Plot to see what the read count/16S copy count relationship looks like.
plot(allqpcr.1$count, allqpcr.1$copies_16S)
#Ok cool, no real relationship here. That's nice, so just seems like fewer



#Subset for only samples with more than 3000 reads. This is what's in the rarefied data.
qpcr3000 <- allqpcr.1[allqpcr.1$count>3000,]
nrow(allqpcr.1)

qpcr3000$div_day <- factor(qpcr3000$div_day, levels=c("LO_D-5", "LO_D1",  "LO_D10","LO_D14" ,"LO_D28", "LO_D31", "LO_D41", 
                                                      "HI_D-5", "HI_D1", "HI_D10","HI_D14" ,"HI_D28", "HI_D31" ,"HI_D41"))

qpcr3000$day_num <-as.numeric(gsub("D", "", qpcr3000$day))
qpcr3000$microb <- factor(qpcr3000$microb, levels=c("LO", "HI"))


# ##### Subset the dfs for subgroups of interest. #####

#Subset to take away all of the samples from the 18684 run. 
qpcr.exp <- qpcr3000[qpcr3000$sac_day != "non_experiment" & qpcr3000$sample_type=="sample",] #all experimental

qpcr.worm <- qpcr3000[qpcr3000$treatment == "worm" & qpcr3000$sample_type != "presample",] #all worm-treated
qpcr.wormhi <- qpcr.worm[qpcr.worm$microb =="HI",]
qpcr.wormlo <- qpcr.worm[qpcr.worm$microb =="LO",]

qpcr.sham <- qpcr3000[qpcr3000$treatment == "control" & qpcr3000$sample_type != "presample",] #all sham-treated
qpcr.shamhi <- qpcr.sham[qpcr.sham$microb =="HI",]
qpcr.shamlo <- qpcr.sham[qpcr.sham$microb =="LO",]

qpcr.preexp <- subset(qpcr3000, sample_type %in% c("test", "presample"))

qpcr.lo <-  qpcr3000[qpcr3000$microb == "LO"  & 
                       qpcr3000$sample_type == "sample",] #lo experimental samples only

qpcr.hi <- qpcr3000[qpcr3000$microb == "HI" & qpcr3000$sample_type == "sample",] #hi experimental samples only


qpcr.imp <- subset(qpcr3000, day %in% c("D-5", "D14", "D41"))
qpcr.imp$MTD <- paste(qpcr.imp$microb, qpcr.imp$day_num, qpcr.imp$treatment, sep="_")
qpcr.imp$MTD <- factor(qpcr.imp$MTD, levels=c("LO_-5_worm","LO_14_worm", "LO_41_worm",
                                              "LO_-5_control", "LO_14_control", "LO_41_control" , 
                                              "HI_-5_worm", "HI_14_worm", "HI_41_worm",  
                                              "HI_-5_control","HI_14_control", "HI_41_control"))

View(qpcr.imp)

# ##### Plot the outcomes of the subsetted data #####


#   ##### Setting colors #####

#Colors, in order
locols <- c('#FDF8FF' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
hicols <- c("#FFFCF3" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")

allcols <- c('#FDF8FF' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049',
             "#FFFCF3" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , '#644C00')


#Need to make sure the boxplot for D-5 can be seen, so change the first color to black. 
otherlo <-  c('#999999' , '#ECC9FE' , '#C66EF5' , '#A11FE5' , '#7E05BE' , '#5C008D' , '#300049')
otherhi <- c("#999999" , "#FFEFBA" , "#FDD964" , "#DCAE1C" , "#BC8F00" , "#8F6D00" , "#644C00")

#psimp timepoint colors
loimpcols <- c('#FDF8FF' , '#A11FE5', '#300049')
hiimpcols <- c("#FFFCF3" , "#DCAE1C" , "#644C00")
impcols <- c('#FDF8FF' , '#A11FE5', '#300049', "#FFFCF3" , "#DCAE1C" , "#644C00")
MTDcols <- c('#FDF8FF' , '#A11FE5' , '#300049', #For use when defining "fill" with 3 timepoints of treatment and controls. 
             '#FFFFFF', '#FFFFFF', '#FFFFFF', 
             "#FFFCF3" , "#DCAE1C" , "#644C00", 
             '#FFFFFF', '#FFFFFF', '#FFFFFF')
MTDcols1 <- c('#ECC9FE' , '#A11FE5' , '#300049', #For use when defining colors with 3 timepoints of treatment and controls. 
              '#ECC9FE' , '#A11FE5' , '#300049',
              "#FFEFBA" , "#DCAE1C" , "#644C00", 
              "#FFEFBA" , "#DCAE1C" , "#644C00" )


precols <- c('#FDF8FF', "#FFFCF3")
d14cols <- c('#A11FE5', "#DCAE1C")
d41cols <- c('#300049', "#644C00")

#   ##### Plotting 16S copies at each day. #####


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



grid.arrange(shamlo, shamhi, ncol=2)


#      ##### Basic stats on the qpcr outcomes. ##### 
mean(qpcr.shamlo$copies_16S)

sd(qpcr.wormlo$copies_16S)
sd(qpcr.wormhi$copies_16S)



#   ##### Plotting 16S copies / number of reads, standardizing for sample quality? #####

qpcr.exp$div_day <- factor(qpcr.exp$div_day, levels=c("LO_D-5", "LO_D1",  "LO_D10","LO_D14" ,"LO_D28", "LO_D31", "LO_D41", 
                                                      "HI_D-5", "HI_D1", "HI_D10","HI_D14" ,"HI_D28", "HI_D31" ,"HI_D41"))
exp <- ggplot(qpcr.exp, aes(x=day, y=copies_16S_per_ng)) +
  geom_boxplot(aes(color=div_day), outlier.shape=NA) + 
  geom_point(aes(shape=treatment, fill=div_day), position=position_jitterdodge(jitter.width=1), alpha=0.8, size=6) + 
  ggtitle("16S copies / DNA conc, all experimental") +
  theme_bw() +
  theme(legend.position="none", text=element_text(size=24)) +
  facet_wrap(~treatment)
exp

#


#Sham- and worm-treated show the same pattern of 16S_copies. Dip in 16S copies and 
# around Day10, but is variable by HI or LO. Observed ASVs also shows a similar pattern. 
# Seems like the samples for Day10 or something are just lower quality. Ask about this in
# the Round lab meeting. 

grid.arrange(shamlo, shamhi, ncol=2)



##### Looking at read count and 16S copies #####
ggplot(qpcr.worm, aes(x=count, y=copies_16S)) +
  geom_point() +
  facet_wrap(~day)


#Just look at the read count. Does it mirror the pattern seen in copies with low 
# points at Day 10-14?
qpcr.worm$div_day <- factor(qpcr.worm$div_day, levels=c("LO_D-5", "LO_D1",  "LO_D10","LO_D14" ,"LO_D28", "LO_D31", "LO_D41", 
                                                        "HI_D-5", "HI_D1", "HI_D10","HI_D14" ,"HI_D28", "HI_D31" ,"HI_D41"))
worm <- ggplot(qpcr.worm, aes(x=day, y=count)) +
  geom_boxplot(aes(color=div_day), outlier.shape=NA) + 
  geom_point(aes(shape=treatment, fill=div_day), position=position_jitterdodge(jitter.width=1), alpha=0.8, size=6) + 
  ggtitle("16S reads from Miseq by day, all worm-treated") +
  scale_fill_manual(values=allcols) +
  scale_shape_manual(values=21) +
  scale_color_manual(values=c("#000000", "#000000","#000000","#000000","#000000","#000000",
                              "#000000","#000000","#000000","#000000","#000000","#000000",
                              "#000000","#000000")) +
  theme_bw() +
  theme(legend.position="none", text=element_text(size=24)) +
  ylim(0, 100000)
worm

#


##### Stats on psimp timepoints #####


library(nlme)
mem.load <- lme(copies_16S ~ microb + treatment + day, random = ~1|mouse, data=qpcr.imp)
anova(mem.load)

emmeans(mem.load, pairwise ~ day)


qpcr.implo <- qpcr.imp[qpcr.imp$microb=="LO",]
qpcr.imphi <- qpcr.imp[qpcr.imp$microb=="HI",]



#
#   ##### Stats on sham-treated during each psimp timepoint #####

qpcrpre <- qpcr.imp[qpcr.imp$day=="D-5" & qpcr.imp$treatment=="control",]
mod <- lm(copies_16S ~ microb + sex, qpcrpre)
summary(mod)


qpcr14 <- qpcr.imp[qpcr.imp$day=="D14" & qpcr.imp$treatment=="control",]
mod <- lm(copies_16S ~ microb + sex, qpcr14) #Outcome is the same with glm. 
summary(mod)


#
##### Create a model for worm-treated  ##### 
#But not sure this is relevant since it's clear that there are patterns across the 
# days that aren't due to parasitism. 
hist(qpcr.worm$copiesbyreads)
hist(qpcr.exp$copies_16S)


qpcr.worm$day_num <- as.numeric(gsub("D", "", qpcr.worm$day))

mod3 <- glm(copies_16S ~ microb + day_num + sex + microb*day_num, qpcr.worm, family=Gamma)
mod.0 <- glm(copies_16S ~ 1, qpcr.worm, family=Gamma)
summary(mod3)
anova(mod3, mod.0, test="F")


#How about modeling copiesbyreads?
mod <- glm(copiesbyreads ~ microb + day_num + sex, qpcr.worm, family=Gamma)
mod1 <- glm(copiesbyreads ~ microb + day_num, qpcr.worm, family=Gamma)
mod2 <- glm(copiesbyreads ~ microb + day_num + microb*day_num, qpcr.worm, family=Gamma)
mod3 <- glm(copiesbyreads ~ microb + day_num + sex + microb*day_num, qpcr.worm, family=Gamma)
anova(mod, mod1, mod2, mod3) #Not much difference at all here either.
AIC(mod2) #AIC values of mod1 is infinite, so shouldn't use it
# than mod2
#



##### Create a model including both worm and sham-treated. #####
hist(qpcr3000$copies_16S)
qpcr3000$day_num <- as.numeric(gsub("D", "", qpcr3000$day))

qpcr.exp$day_num <- as.numeric(gsub("D", "", qpcr.exp$day))

mod <- glm(copies_16S ~ treatment + microb + day_num + sex + microb*day_num + treatment*day_num, data = qpcr.exp, family = Gamma)
mod.0 <- glm(copies_16S ~ 1, qpcr.exp, family=Gamma)
anova(mod, mod.0, test = "F")
summary(mod)



#




