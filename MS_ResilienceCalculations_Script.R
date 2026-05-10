#Pulling together resilience calculations. 

#Creating dfs for the calculation. 
#Load workspace with the correct objects. 
load("/Users/mdoolin/Documents/1U4U_16S/GG2_reclassification/23Feb2026_Whole1U4U_16SWorkspace.RData")

#Bring in metadata for psexp3000no31. Already pulled it in alphabetachange script. 
dfmetadat 

###### Test whether worm-treated group changed between 28 dpi and 41 dpi. #####
# (This part is the same approach we used to test group against itself between other timepoints.)

#   ##### 1. Subset the data to each group. #####

#Dfs with alpha diversity
dflotreat <- dfmetadat[dfmetadat$MT=="LO_worm",]
dfhitreat <- dfmetadat[dfmetadat$MT=="HI_worm",]
#double check that the correct rows were pulled.
unique(dflotreat$mouse)

#ps objects for beta diversity calculations
library(phyloseq)
pslotreat <- subset_samples(psexp3000no31, microb == "LO" & treatment == "worm")
pslotreat <- prune_taxa(taxa_sums(pslotreat) >= 1, pslotreat)

pshitreat <- subset_samples(psexp3000no31, microb == "HI" & treatment == "worm")
pshitreat <- prune_taxa(taxa_sums(pshitreat) >= 1, pshitreat)


#   ##### 2. Further subset to two timepoints of interest.#####
#      ##### df's for alpha diversity #####
dflotreatpre1 <- subset(dflotreat, day %in% c('D-5', 'D1')) #64 rows 
dfhitreatpre1 <- subset(dfhitreat, day %in% c('D-5', 'D1')) #36 rows 

dflotreatpre10 <- subset(dflotreat, day %in% c('D-5', 'D10'))
dfhitreatpre10 <- subset(dfhitreat, day %in% c('D-5', 'D10'))

dflotreatpre14 <- subset(dflotreat, day %in% c('D-5', 'D14'))
dfhitreatpre14 <- subset(dfhitreat, day %in% c('D-5', 'D14'))

dflotreatpre28 <- subset(dflotreat, day %in% c('D-5', 'D28'))
dfhitreatpre28 <- subset(dfhitreat, day %in% c('D-5', 'D28'))

dflotreatpre41 <- subset(dflotreat, day %in% c('D-5', 'D41'))
dfhitreatpre41 <- subset(dfhitreat, day %in% c('D-5', 'D41'))

#Starting at timepoints other than pretreatment. 
dflotreat110 <- subset(dflotreat, day %in% c('D1', 'D10'))
dfhitreat110 <- subset(dfhitreat, day %in% c('D1', 'D10'))

dflotreat1014 <- subset(dflotreat, day %in% c('D14', 'D10'))
dfhitreat1014 <- subset(dfhitreat, day %in% c('D14', 'D10'))

dflotreat1428 <- subset(dflotreat, day %in% c('D14', 'D28'))
dfhitreat1428 <- subset(dfhitreat, day %in% c('D14', 'D28'))

dflotreat1441 <- subset(dflotreat, day %in% c('D14', 'D41'))
dfhitreat1441 <- subset(dfhitreat, day %in% c('D14', 'D41'))

dflotreat2841 <- subset(dflotreat, day %in% c('D28', 'D41'))
#38 rows, and not all paired. 
dfhitreat2841 <- subset(dfhitreat, day %in% c('D28', 'D41'))
#Only 13 rows, and not all paired. 


#      ##### ps objects for beta diversity #####
# first for LO
pslotreatpre1 <- subset_samples(pslotreat, day %in% c("D-5", "D1"))
pslotreatpre1 <- prune_taxa(taxa_sums(pslotreatpre1) >= 1, pslotreatpre1)

pslotreatpre10 <- subset_samples(pslotreat, day %in% c("D-5", "D10"))
pslotreatpre10 <- prune_taxa(taxa_sums(pslotreatpre10) >= 1, pslotreatpre10)

pslotreatpre14 <- subset_samples(pslotreat, day %in% c("D-5", "D14"))
pslotreatpre14 <- prune_taxa(taxa_sums(pslotreatpre14) >= 1, pslotreatpre14)

pslotreatpre28 <- subset_samples(pslotreat, day %in% c("D-5", "D28"))
pslotreatpre28 <- prune_taxa(taxa_sums(pslotreatpre28) >= 1, pslotreatpre28)

pslotreatpre41 <- subset_samples(pslotreat, day %in% c("D-5", "D41"))
pslotreatpre41 <- prune_taxa(taxa_sums(pslotreatpre41) >= 1, pslotreatpre41)

pslotreat110 <- subset_samples(pslotreat, day %in% c("D1", "D10"))
pslotreat110 <- prune_taxa(taxa_sums(pslotreat110) >= 1, pslotreat110)

pslotreat1014 <- subset_samples(pslotreat, day %in% c("D14", "D10"))
pslotreat1014 <- prune_taxa(taxa_sums(pslotreat1014) >= 1, pslotreat1014)

pslotreat1428 <- subset_samples(pslotreat, day %in% c("D14", "D28"))
pslotreat1428 <- prune_taxa(taxa_sums(pslotreat1428) >= 1, pslotreat1428)

pslotreat1441 <- subset_samples(pslotreat, day %in% c("D14", "D41"))
pslotreat1441 <- prune_taxa(taxa_sums(pslotreat1441) >= 1, pslotreat1441)

pslotreat2841 <- subset_samples(pslotreat, day %in% c("D28", "D41"))
pslotreat2841 <- prune_taxa(taxa_sums(pslotreat2841) >= 1, pslotreat2841)


#Then for HI. 
pshitreatpre1 <- subset_samples(pshitreat, day %in% c("D-5", "D1"))
pshitreatpre1 <- prune_taxa(taxa_sums(pshitreatpre1) >= 1, pshitreatpre1)

pshitreatpre10 <- subset_samples(pshitreat, day %in% c("D-5", "D10"))
pshitreatpre10 <- prune_taxa(taxa_sums(pshitreatpre10) >= 1, pshitreatpre10)

pshitreatpre14 <- subset_samples(pshitreat, day %in% c("D-5", "D14"))
pshitreatpre14 <- prune_taxa(taxa_sums(pshitreatpre14) >= 1, pshitreatpre14)

pshitreatpre28 <- subset_samples(pshitreat, day %in% c("D-5", "D28"))
pshitreatpre28 <- prune_taxa(taxa_sums(pshitreatpre28) >= 1, pshitreatpre28)

pshitreatpre41 <- subset_samples(pshitreat, day %in% c("D-5", "D41"))
pshitreatpre41 <- prune_taxa(taxa_sums(pshitreatpre41) >= 1, pshitreatpre41)

pshitreat110 <- subset_samples(pshitreat, day %in% c("D1", "D10"))
pshitreat110 <- prune_taxa(taxa_sums(pshitreat110) >= 1, pshitreat110)

pshitreat1014 <- subset_samples(pshitreat, day %in% c("D14", "D10"))
pshitreat1014 <- prune_taxa(taxa_sums(pshitreat1014) >= 1, pshitreat1014)

pshitreat1428 <- subset_samples(pshitreat, day %in% c("D14", "D28"))
pshitreat1428 <- prune_taxa(taxa_sums(pshitreat1428) >= 1, pshitreat1428)

pshitreat1441 <- subset_samples(pshitreat, day %in% c("D14", "D41"))
pshitreat1441 <- prune_taxa(taxa_sums(pshitreat1441) >= 1, pshitreat1441)

pshitreat2841 <- subset_samples(pshitreat, day %in% c("D28", "D41"))
pshitreat2841 <- prune_taxa(taxa_sums(pshitreat2841) >= 1, pshitreat2841)
#
#

#   ##### 3. Create the model looking for change. #####
#Can use a poisson if mean = variance. But use negative binomial if mean not equal to variance. 
# Function to check for overdispersion (can be sourced from online resources like the QCBS R workshop wiki)
# visit the tutoial for more exercises: https://r.qcbs.ca/workshop07/book-en/poisson-glmm.html
# Install and load necessary packages
# install.packages("lme4")
library(lme4)
# install.packages("glmmTMB")
library(glmmTMB)



#      ##### Richness #####


#         ##### pretreatment to D1 #####
# For the LO group #
# Fit the Poisson GLMM. Key assumption is that the mean == variance. 
mean(dflotreatpre1$Observed) #30.32812
var(dflotreatpre1$Observed) #144.097
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreatpre1, 
                 family = nbinom1) # or nbinom2
# nbinom2 is for when variance grows faster than mean. nbonim1 is when variance grows linearly with mean. 
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  3.31377    0.06225   53.23   <2e-16 ***
##dayD1        0.09797    0.04545    2.16   0.0311 *  

# For the HI group #
mean(dfhitreatpre1$Observed) #166.3611
var(dfhitreatpre1$Observed) #1001.552
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.08735    0.04384  116.05   <2e-16 ***
##dayD1        0.04767    0.05592    0.85    0.394 

#

ggplot(dflotreatpre1, aes(x=day, y=Observed)) + 
  geom_boxplot() +
  geom_point(aes(color=mouse)) +
  theme_bw()

#
#         ##### pretreatment to D10 #####
# For the LO group #
mean(dflotreatpre10$Observed) #46.76562
var(dflotreatpre10$Observed) #479.9601
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreatpre10, 
                 family = nbinom1) 
summary(lomod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  3.35084    0.05759   58.18   <2e-16 ***
##dayD10       0.81799    0.06618   12.36   <2e-16 ***

# For the HI group #
mean(dfhitreatpre10$Observed) #162.5833
var(dfhitreatpre10$Observed) #1263.793
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreatpre10, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.095762   0.050907  100.10   <2e-16 ***
##dayD10      -0.009164   0.071835   -0.13    0.898   

#
#         ##### pretreatment to D14 #####
# For the LO group #
mean(dflotreatpre14$Observed) #44.51562
var(dflotreatpre14$Observed) #367.9045
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreatpre14, 
                 family = nbinom1) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  4.16244    0.04690   88.75   <2e-16 ***
##dayD41      -0.05720    0.04477   -1.28    0.201    

# For the HI group #
mean(dfhitreatpre14$Observed) #163.5152
var(dfhitreatpre14$Observed) #971.0076
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreatpre14, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.09024    0.04371  116.46   <2e-16 ***
##dayD14       0.01460    0.06437    0.23    0.821   

#
#         ##### pretreatment to D28 #####
# For the LO group #
mean(dflotreatpre28$Observed) #39.93617
var(dflotreatpre28$Observed) #386.148
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreatpre28, 
                 family = nbinom1) 
summary(lomod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  3.35190    0.05687   58.93   <2e-16 ***
##dayD28       0.81029    0.07907   10.25   <2e-16 ***

# For the HI group #
mean(dfhitreatpre28$Observed) #148.6957
var(dfhitreatpre28$Observed) #1527.949
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreatpre28, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.08947    0.04257   119.6  < 2e-16 ***
##dayD28      -0.48723    0.11073    -4.4 1.08e-05 ***

#
#         ##### pretreatment to D41 #####
# For the LO group #
mean(dflotreatpre41$Observed) #42.41818
var(dflotreatpre41$Observed) #417.2108
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreatpre41, 
                 family = nbinom1) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  3.32046    0.05986   55.47   <2e-16 ***
##dayD41       0.80804    0.04664   17.32   <2e-16 *** 

# For the HI group #
mean(dfhitreatpre41$Observed) #159.1923
var(dfhitreatpre41$Observed) #1019.682
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreatpre41, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.08028    0.04451  114.13   <2e-16 ***
##dayD41      -0.07357    0.06053   -1.22    0.224 

#

#         ##### D1 to D10 #####
# For the LO group #
hist(dflotreat110$Observed) #48.25
mean(dflotreat110$Observed) #48.25
var(dflotreat110$Observed) #430.9524
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreat110, 
                 family = nbinom1) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  3.45068    0.05454   63.27   <2e-16 ***
##dayD10       0.71836    0.06375   11.27   <2e-16 ***

# For the HI group #
mean(dfhitreat110$Observed) #166.6111
var(dfhitreat110$Observed) #1363.159
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreat110, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.13879    0.05116   100.5   <2e-16 ***
##dayD10      -0.05349    0.06715    -0.8    0.426    

#


#         ##### D10 to D14 #####
# For the LO group #
hist(dflotreat1014$Observed)
mean(dflotreat1014$Observed) #62.4375
var(dflotreat1014$Observed) #138.1548
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreat1014, 
                 family = nbinom2) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  4.15942    0.03344  124.40   <2e-16 ***
##dayD14      -0.07210    0.03166   -2.28   0.0228 *  

# For the HI group #
mean(dfhitreat1014$Observed) #163.7879
var(dfhitreat1014$Observed) #1368.11
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreat1014, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.08651    0.05272   96.48   <2e-16 ***
##dayD14       0.02130    0.07462    0.29    0.775   

#



#         ##### D14 to D28 #####
# For the LO group #
mean(dflotreat1428$Observed) #61.2766
var(dflotreat1428$Observed) #92.24792
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreat1428, 
                 family = nbinom1) 
summary(lomod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  4.09262    0.02875   142.4   <2e-16 ***
##dayD28       0.05506    0.04224     1.3    0.192   

# For the HI group #
mean(dfhitreat1428$Observed) #148.6
var(dfhitreat1428$Observed) #1775.516
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreat1428, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.09620    0.04948  102.99  < 2e-16 ***
##dayD28      -0.48774    0.09844   -4.95 7.24e-07 ***

#
#         ##### D14 to D41 #####
# For the LO group #
mean(dflotreat1441$Observed) #60.65455
var(dflotreat1441$Observed) #132.3044
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreat1441, 
                 family = nbinom1) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  4.09304    0.03347  122.28   <2e-16 ***
##dayD41       0.02482    0.04442    0.56    0.576     

# For the HI group #
mean(dfhitreat1441$Observed) #160.4783
var(dfhitreat1441$Observed) #1149.261
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreat1441, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  5.09177    0.05105   99.75   <2e-16 ***
##dayD41      -0.05848    0.07114   -0.82    0.411   

#
#         ##### D28 to D41 #####
# For the LO group #
mean(dflotreat2841$Observed) #62.21053
var(dflotreat2841$Observed) #138.9815
lomod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dflotreat2841, 
                 family = nbinom1) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  4.16244    0.04690   88.75   <2e-16 ***
##dayD41      -0.05720    0.04477   -1.28    0.201    

# For the HI group #
mean(dfhitreat2841$Observed) #131.9231
var(dfhitreat2841$Observed) #1676.244
himod <- glmmTMB(Observed ~ day + (1 | mouse), 
                 data = dfhitreat2841, 
                 family = nbinom1) # or nbinom1. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  4.69268    0.10482   44.77  < 2e-16 ***
##dayD41       0.31453    0.09676    3.25  0.00115 ** 

#
#      ##### Evenness #####
#         ##### pretreatment to D1 #####
# For the LO group #
# Fit the Poisson GLMM. Key assumption is that the mean == variance. 
mean(dflotreatpre1$evenness) #0.5792791
var(dflotreatpre1$evenness) #0.003284214
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreatpre1, 
                 family = gaussian) 
# nbinom2 is for when variance grows faster than mean. nbonim1 is when variance grows linearly with mean. 
# View model summary
summary(lomod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.587926   0.009934   59.18   <2e-16 ***
##dayD1       -0.017294   0.010892   -1.59    0.112    

# For the HI group #
mean(dfhitreatpre1$evenness) #0.7245048
var(dfhitreatpre1$evenness) #0.001006412
# Fit the negative binomial GLMM
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreatpre1, 
                 family = gaussian) # or gaussian. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.723544   0.007369   98.18   <2e-16 ***
##dayD1       0.001921   0.010422    0.18    0.854   

#
#         ##### pretreatmenåt to D10 #####
# For the LO group #
mean(dflotreatpre10$evenness) #0.5908121
var(dflotreatpre10$evenness) #0.00583631
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreatpre10, 
                 family = gaussian) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.587926   0.013389   43.91   <2e-16 ***
##dayD10      0.005772   0.017191    0.34    0.737 

# For the HI group #
mean(dfhitreatpre10$evenness) #0.7142809
var(dfhitreatpre10$evenness) #0.001455466
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreatpre10, 
                 family = gaussian) # or gaussian. 
summary(himod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.723544   0.008593   84.20   <2e-16 ***
##dayD10      -0.018527   0.012153   -1.52    0.127  

#
#         ##### pretreatment to D14 #####
# For the LO group #
mean(dflotreatpre14$evenness) #0.579983
var(dflotreatpre14$evenness) #0.005278023
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreatpre14, 
                 family = gaussian) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.58793    0.01266   46.42   <2e-16 ***
##dayD14      -0.01589    0.01791   -0.89    0.375   

# For the HI group #
mean(dfhitreatpre14$evenness) #0.7159633
var(dfhitreatpre14$evenness) #0.001162381
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreatpre14, 
                 family = gaussian) # or gaussian. 
summary(himod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.723544   0.007667   94.38   <2e-16 ***
##dayD14      -0.016553   0.010920   -1.52     0.13   

#
#         ##### pretreatment to D28 #####
# For the LO group #
mean(dflotreatpre28$evenness) #0.585516
var(dflotreatpre28$evenness) #0.004906898
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreatpre28, 
                 family = gaussian) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.587926   0.012235   48.05   <2e-16 ***
##dayD28      -0.007399   0.021593   -0.34    0.732   

# For the HI group #
mean(dfhitreatpre28$evenness) #0.7259308
var(dfhitreatpre28$evenness) #0.0005588508
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreatpre28, 
                 family = gaussian) # or gaussian. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.723544   0.005344  135.39   <2e-16 ***
##dayD28      0.010978   0.011462    0.96    0.338  

#
#         ##### pretreatment to D41 #####
# For the LO group #
mean(dflotreatpre41$evenness) #0.5781374
var(dflotreatpre41$evenness) #0.00453715
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreatpre41, 
                 family = gaussian) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.58793    0.01161   50.65   <2e-16 ***
##dayD41      -0.02360    0.01564   -1.51    0.131  

# For the HI group #
mean(dfhitreatpre41$evenness) #0.7243849
var(dfhitreatpre41$evenness) #0.0005633282
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreatpre41, 
                 family = gaussian) # or gaussian. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.723544   0.005477  132.11   <2e-16 ***
##dayD41      0.002590   0.009744    0.27     0.79

#
#         ##### D1 to D10 #####
# For the LO group #
mean(dflotreat110$evenness) #0.5821649
sd(dflotreat110$evenness) #0.07774749
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreat110, 
                 family = gaussian) 
summary(lomod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.57063    0.01348   42.32   <2e-16 ***
##dayD10       0.02307    0.01889    1.22    0.222     


# For the HI group #
mean(dfhitreat110$evenness) #0.7152414
sd(dfhitreat110$evenness) #0.04310924
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreat110, 
                 family = gaussian) # or nbinom1. 
summary(himod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.725465   0.009725   74.60   <2e-16 ***
##dayD10      -0.020448   0.011027   -1.85   0.0637 .

#

#         ##### D10 to D14 #####
# For the LO group #
mean(dflotreat1014$evenness) #0.5828688
sd(dflotreat1014$evenness) #0.08963456
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreat1014, 
                 family = gaussian) 
summary(lomod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.59370    0.01560   38.05   <2e-16 ***
##dayD14      -0.02166    0.02148   -1.01    0.313  

# For the HI group #
mean(dfhitreat1014$evenness) #0.7058578
sd(dfhitreat1014$evenness) #0.04367432
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreat1014, 
                 family = gaussian) # or nbinom1. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.705018   0.010135   69.57   <2e-16 ***
##dayD14      0.001848   0.015032    0.12    0.902  

#
#         ##### D14 to D28 #####
# For the LO group #
mean(dflotreat1428$evenness) #0.5746996
var(dflotreat1428$evenness) #0.00776845
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreat1428, 
                 family = gaussian) 
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.572040   0.015399   37.15   <2e-16 ***
##dayD28      0.008334   0.027258    0.31     0.76   

# For the HI group #
mean(dfhitreat1428$evenness) #0.71378
var(dfhitreat1428$evenness) #0.001465767
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreat1428, 
                 family = gaussian) # or gaussian. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.706866   0.009125   77.46   <2e-16 ***
##dayD28      0.027656   0.018251    1.52     0.13    

#
#         ##### D14 to D41 #####
# For the LO group #
mean(dflotreat1441$evenness) #0.5688943
var(dflotreat1441$evenness) #0.006850653
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreat1441, 
                 family = gaussian) 
summary(lomod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.572040   0.014483   39.50   <2e-16 ***
##dayD41      -0.007522   0.022396   -0.34    0.737  

# For the HI group #
mean(dfhitreat1441$evenness) #0.7136174
var(dfhitreat1441$evenness) #0.001325717
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreat1441, 
                 family = gaussian) # or gaussian. 
summary(himod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 0.706325   0.008716   81.03   <2e-16 ***
##dayD41      0.024367   0.011226    2.17     0.03 *  

#

#         ##### D28 to D41 #####

# For the LO group #
# Fit the Poisson GLMM. Key assumption is that the mean == variance. 
mean(dflotreat2841$evenness) #0.5707771
var(dflotreat2841$evenness) #0.007246438
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dflotreat2841, 
                 family = gaussian) # or nbinom2
# nbinom2 is for when variance grows faster than mean. nbonim1 is when variance grows linearly with mean. 
# View model summary
summary(lomod)
##           Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.58886    0.02384  24.697   <2e-16 ***
##dayD41      -0.02434    0.02710  -0.898    0.369    

# For the HI group #
mean(dfhitreat2841$evenness) #0.7294478
var(dfhitreat2841$evenness) #0.0003593444
# Fit the negative binomial GLMM
himod <- glmmTMB(evenness ~ day + (1 | mouse), 
                 data = dfhitreat2841, 
                 family = gaussian) # or gaussian. 
summary(himod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  0.734522   0.007945   92.45   <2e-16 ***
##dayD41      -0.008246   0.010128   -0.81    0.416  




#      ##### qPCR #####


#         ##### First, a little qc on the qPCR data to compare HI and LO loads at pretreatment #####
colnames(dfmetadat)

#Look at the relationship between the qPCR input DNA, copies of 16S, and day
boxplot(dfmetadat$qPCR_input_DNA ~ dfmetadat$day, main = "qPCR input DNA, per day")
boxplot(dfmetadat$copies_16S_per_ng ~ dfmetadat$day, main = "Copies 16S per ng DNA, per day")

#Before looking at change, just look at the raw data at pretreatment, D14, and D41 to see if there is a difference among groups. 
colnames(dfmetadat)
colnames(dfpre)
hist(dfpre$copies_16S_per_ng.pre)
dftmppre <- dfpre %>% drop_na(copies_16S_per_ng.pre)
boxplot(dftmppre$copies_16S_per_ng.pre ~ dftmppre$MT)
summary(dftmppre$copies_16S_per_ng.pre) 
##    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 2211564  3817133  4254122  4406665  4759495 13264163 
sd(dftmppre$copies_16S_per_ng.pre) #1489552
premod <- glm.nb(copies_16S_per_ng.pre ~ microb, data = dftmppre) 
# nbinom2 is for when variance grows faster than mean. nbonim1 is when variance grows linearly with mean. 
# View model summary
summary(premod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.34581    0.04375 350.750   <2e-16 ***
##microbHI    -0.14599    0.07511  -1.944   0.0519 .  
#Borderline. This is actually a result of there being more reads in LO animals. 
#
lomod <- glm.nb(copies_16S_per_ng.pre ~ treatment, data = dftmppre[dftmppre$microb=="LO",]) 
summary(lomod)
##              Estimate Std. Error z value Pr(>|z|)    
##(Intercept)    15.2418     0.1178 129.374   <2e-16 ***
##treatmentworm   0.1229     0.1287   0.955     0.34   

himod <- glm.nb(copies_16S_per_ng.pre ~ treatment, data = dftmppre[dftmppre$microb=="HI",]) 
summary(himod)
##              Estimate Std. Error z value Pr(>|z|)    
##(Intercept)    14.7521     0.1072  137.67  < 2e-16 ***
##treatmentworm   0.4893     0.1133    4.32 1.56e-05 ***
#This significant difference is a result of the HI sham having low loads pretreatment. Weird.  a

#         ##### pretreatment to D1 #####
# For the LO group #
# Fit the GLMM. Key assumption is that the mean == variance if using Poisson, so we won't use that. 
dftmplo <- dflotreatpre1 %>% drop_na(copies_16S_per_ng)

mean(dftmplo$copies_16S_per_ng) #4279216
var(dftmplo$copies_16S_per_ng) #2.172715e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# nbinom2 is for when variance grows faster than mean. nbonim1 is when variance grows linearly with mean. 
# View model summary
summary(lomod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.35178    0.04851   316.5  < 2e-16 ***
##dayD1       -0.18499    0.06179    -3.0  0.00275 ** 
# For the HI group #
dftmphi <- dfhitreatpre1 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #3034234
var(dftmphi$copies_16S_per_ng) #11.810707e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.24141    0.07246  210.34  < 2e-16 ***
##dayD1       -0.78042    0.10248   -7.62 2.62e-14 ***

#
#         ##### pretreatment to D10 #####
# For the LO group #
# Fit the GLMM. Key assumption is that the mean == variance if using Poisson, so we won't use that. 
dftmplo <- dflotreatpre10 %>% drop_na(copies_16S_per_ng)

mean(dftmplo$copies_16S_per_ng) #3317944
var(dftmplo$copies_16S_per_ng) #3.780618e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.36475    0.05628  273.00   <2e-16 ***
##dayD10      -0.87047    0.07897  -11.02   <2e-16 ***

# For the HI group #
dftmphi <- dfhitreatpre10 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #3034234
var(dftmphi$copies_16S_per_ng) #11.810707e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  15.2418     0.1021  149.33   <2e-16 ***
##dayD10       -0.3387     0.1520   -2.23   0.0259 *  

#
#         ##### pretreatment to D14 #####
# For the LO group #
# Fit the GLMM. Key assumption is that the mean == variance if using Poisson, so we won't use that. 
dftmplo <- dflotreatpre14 %>% drop_na(copies_16S_per_ng)

mean(dftmplo$copies_16S_per_ng) #3659888
var(dftmplo$copies_16S_per_ng) #6.936947e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.38407    0.09191  167.38  < 2e-16 ***
##dayD14      -0.71674    0.11105   -6.45 1.09e-10 ***

# For the HI group #
dftmphi <- dfhitreatpre14 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #3594537
var(dftmphi$copies_16S_per_ng) #1.048456e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.24166    0.06129  248.66  < 2e-16 ***
##dayD14      -0.36525    0.09339   -3.91 9.19e-05 ***

#
#         ##### pretreatment to D28 #####
# For the LO group #
# Fit the GLMM. Key assumption is that the mean == variance if using Poisson, so we won't use that. 
dftmplo <- dflotreatpre28 %>% drop_na(copies_16S_per_ng)
#Remove the lo outlier, which is definitely a problem with only 7 copies per ng. 
dftmplo <- dftmplo[dftmplo$copies_16S_per_ng > 100, ]

mean(dftmplo$copies_16S_per_ng) #3987477
var(dftmplo$copies_16S_per_ng) #3.650667e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.34871    0.05372  285.72  < 2e-16 ***
##dayD28      -0.65290    0.08742   -7.47 8.12e-14 ***
#

# For the HI group #
dftmphi <- dfhitreatpre28 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #4116502
var(dftmphi$copies_16S_per_ng) #508256056330
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.23659    0.04309   353.6   <2e-16 ***
##dayD28      -0.07901    0.09475    -0.8    0.404    

#
#         ##### pretreatment to D41 #####
# For the LO group #
dftmplo <- dflotreatpre41 %>% drop_na(copies_16S_per_ng)
#Remove the lo outlier, which is definitely a problem with only 7 copies per ng. 

mean(dftmplo$copies_16S_per_ng) #4780714
var(dftmplo$copies_16S_per_ng) #2.69639e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.34899    0.05305  289.34   <2e-16 ***
##dayD41       0.05425    0.07400    0.73    0.463    
#

# For the HI group #
dftmphi <- dfhitreatpre41 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #4203117
var(dftmphi$copies_16S_per_ng) #441073455741
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.23332    0.03971   383.6   <2e-16 ***
##dayD41       0.04480    0.05085     0.9    0.378    

#
#         ##### D1 to D10 #####
# For the LO group #
# Fit the GLMM. Key assumption is that the mean == variance if using Poisson, so we won't use that. 
dftmplo <- dflotreat110 %>% drop_na(copies_16S_per_ng)

mean(dftmplo$copies_16S_per_ng) #2917700
var(dftmplo$copies_16S_per_ng) #1.453128e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.16635    0.04776   317.6   <2e-16 ***
##dayD10      -0.67925    0.06434   -10.6   <2e-16 *** 
#

# For the HI group #
dftmphi <- dfhitreat110 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #2446194
var(dftmphi$copies_16S_per_ng) #2.085839e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  14.4610     0.1184  122.14  < 2e-16 ***
##dayD10        0.4483     0.1674    2.68  0.00742 ** 

#
#         ##### D10 to D14 #####
# For the LO group #
# Fit the GLMM. Key assumption is that the mean == variance if using Poisson, so we won't use that. 
dftmplo <- dflotreat1014 %>% drop_na(copies_16S_per_ng)

mean(dftmplo$copies_16S_per_ng) #2917700
var(dftmplo$copies_16S_per_ng) #1.453128e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 14.49648    0.08795  164.82   <2e-16 ***
##dayD14       0.14770    0.09646    1.53    0.126 
#

# For the HI group #
dftmphi <- dfhitreat1014 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #2949590
var(dftmphi$copies_16S_per_ng) #2.059999e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 14.90928    0.11581  128.74   <2e-16 ***
##dayD14      -0.02699    0.17233   -0.16    0.876 

#

#         ##### D14 to D28 #####
# For the LO group #
dftmplo <- dflotreat1428 %>% drop_na(copies_16S_per_ng)
#Remove the lo outlier, which is definitely a problem with only 7 copies per ng. 
dftmplo <- dftmplo[dftmplo$copies_16S_per_ng > 100, ]

mean(dftmplo$copies_16S_per_ng) #2567852
var(dftmplo$copies_16S_per_ng) #5.930128e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 14.58754    0.09728  149.95   <2e-16 ***
##dayD28       0.13696    0.09163    1.49    0.135    
#

# For the HI group #
dftmphi <- dfhitreat1428 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #3074528
var(dftmphi$copies_16S_per_ng) #1.083843e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept)  14.8823     0.0857  173.65   <2e-16 ***
##dayD28        0.2842     0.2040    1.39    0.164    

#
#         ##### D14 to D41 #####
# For the LO group #
dftmplo <- dflotreat1441 %>% drop_na(copies_16S_per_ng)
#Remove the lo outlier, which is definitely a problem with only 7 copies per ng. 
dftmplo <- dftmplo[dftmplo$copies_16S_per_ng > 100, ]

mean(dftmplo$copies_16S_per_ng) #3533620
var(dftmplo$copies_16S_per_ng) #6.889424e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 14.58340    0.09672   150.8   <2e-16 ***
##dayD41       0.84883    0.07071    12.0   <2e-16 ***
#

# For the HI group #
dftmphi <- dfhitreat1441 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #3330431
var(dftmphi$copies_16S_per_ng) #1.245248e+12
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 14.84470    0.08565  173.32  < 2e-16 ***
##dayD41       0.29609    0.08475    3.49 0.000477 ***

#
#         ##### D28 to D41 #####
# For the LO group #
dftmplo <- dflotreat2841 %>% drop_na(copies_16S_per_ng)
#Remove the lo outlier, which is definitely a problem with only 7 copies per ng. 
dftmplo <- dftmplo[dftmplo$copies_16S_per_ng > 100, ]

mean(dftmplo$copies_16S_per_ng) #3889867
var(dftmplo$copies_16S_per_ng) #2.723873e+12
hist(dftmplo$copies_16S_per_ng)
boxplot(dftmplo$copies_16S_per_ng ~ dftmplo$day)
#Mean and variance are not equal. Should use negative binomial. 
lomod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmplo, 
                 family = nbinom2) # or nbinom2
# View model summary
summary(lomod)
##            Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 14.57865    0.06923  210.58   <2e-16 ***
##dayD41       0.79898    0.05864   13.62   <2e-16 ***
#

# For the HI group #
dftmphi <- dfhitreat2841 %>% drop_na(copies_16S_per_ng)

mean(dftmphi$copies_16S_per_ng) #4167556
var(dftmphi$copies_16S_per_ng) #718965554298
hist(dftmphi$copies_16S_per_ng)
boxplot(dftmphi$copies_16S_per_ng ~ dftmphi$day)
himod <- glmmTMB(copies_16S_per_ng ~ day + (1 | mouse), 
                 data = dftmphi, 
                 family = nbinom2) # or nbinom1. 
summary(himod)
# Fit the negative binomial GLMM
##             Estimate Std. Error z value Pr(>|z|)    
##(Intercept) 15.26026    0.09885  154.38  < 2e-16 ***
##dayD41      -0.07802    0.02614   -2.98  0.00284 ** 

#

#      ##### Generalized UniFrac #####
#         ##### Plotting by day just to see what's up #####
obj <- pshitreat2841
tree <- phyloseq::phy_tree(obj) #Extract the tree. 
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0, 0.5, 1))$unifracs #Call as many weightings as you want. 
d50 <- gunif[, , "d_0.5"]
GU50pcoa <- data.frame(cmdscale(d50, eig=FALSE))
GU50pcoa <- GU50pcoa %>% 
  dplyr::rename("GU50_pc1" = "X1",
                "GU50_pc2" = "X2")
GU50pcoa <- GU50pcoa*(-1)
GU50pcoa$row <- rownames(GU50pcoa)
dftmp <- merge(dfexp3000no31, GU50pcoa, by="row", all.x=FALSE)

pc <- pcoa(d50) #generate the pcoa from your distance matrix so that you can grab the axis values to label the axes. 
GU50plot <- ggplot(aes(x=GU50_pc1, y=GU50_pc2), data=dftmp) +
  geom_point(aes(color=day), size=3) +
  stat_ellipse(aes(color= day), type="t", level=0.95)+
  #scale_color_manual(values=sexcols) +
  ylim(-0.5, 0.5) + xlim(-0.5, 0.5) +
  theme_bw() +
  xlab(paste("GUnif Axis 1 (", label_percent(accuracy=0.01)(pc$values[1,2]), ")", sep="")) + 
  ylab(paste("GUnif Axis 2 (", label_percent(accuracy=0.01)(pc$values[2,2]), ")", sep="")) +
  theme(text=element_text(size=8)) + 
  ggtitle("GUniFrac (50%) PCOA")
GU50plot

#
#         ##### Running the PERMANOVA #####

library(GUniFrac)

#Generating GUniFrac distances. 
obj <- pslotreat1014 #Which ps object? See below for all that I ran. 
tree <- phyloseq::phy_tree(obj) #Extract the tree. 
otu <- t(otu_table(obj)) #Need to transpose for this package.
gunif <- GUniFrac(otu, tree, alpha=c(0.5))$unifracs #Call as many weightings as you want. 
d50 <- gunif[, , "d_0.5"] #This is pulling the distance matrix of your choice. 
metadata <- as(sample_data(obj), "data.frame")  
mod <-  adonis3(d50 ~ day, strata = metadata$mouse, permutations = 5000, data=metadata) 
mod$aov.tab #View your p-values. 
#
#

#         ##### Outcomes for LO #####
pslotreatpre1
##          Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
##day        1    0.0853 0.085299  2.3435 0.03642 0.0007998 ***
##Residuals 62    2.2567 0.036399         0.96358              
##Total     63    2.3420                  1.00000      

pslotreatpre10
##          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##day        1    1.3754 1.37537   26.93 0.30282  2e-04 ***
##Residuals 62    3.1665 0.05107         0.69718           
##Total     63    4.5418                 1.00000        

pslotreatpre14
##          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##day        1    1.3822  1.3822  27.211 0.30502  2e-04 ***
##Residuals 62    3.1494  0.0508         0.69498           
##Total     63    4.5316                 1.00000     

pslotreatpre28
##         Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##day        1    1.1016 1.10160  27.123 0.37606  2e-04 ***
##Residuals 45    1.8277 0.04062         0.62394           
##Total     46    2.9293                 1.00000    

pslotreatpre41
##          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##day        1    1.1097  1.1097  27.399 0.34079  2e-04 ***
##Residuals 53    2.1466  0.0405         0.65921           
##Total     54    3.2563                 1.00000  

pslotreat110
##          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
##day        1    0.9946 0.99456  17.425 0.21939  2e-04 ***
##Residuals 62    3.5388 0.05708         0.78061           
##Total     63    4.5334                 1.00000   

pslotreat1014
##          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
##day        1    0.1566 0.156646  2.1572 0.03362 0.0124 *
##Residuals 62    4.5023 0.072617         0.96638         
##Total     63    4.6589                  1.00000     

pslotreat1428
##          Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
##day        1    0.1004 0.10042  1.4284 0.03077 0.08098 .
##Residuals 45    3.1635 0.07030         0.96923          
##Total     46    3.2639                 1.00000 

pslotreat1441
##          Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
##day        1    0.1613 0.161314  2.5311 0.04558 0.0009998 ***
##Residuals 53    3.3779 0.063734         0.95442              
##Total     54    3.5392                  1.00000    

pslotreat2841
##          Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)
##day        1   0.08592 0.085918    1.43 0.0382 0.1112
##Residuals 36   2.16297 0.060083         0.9618       
##Total     37   2.24889                  1.0000  


#         ##### Outcomes for HI #####

pshitreatpre1
##          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
##day        1   0.04185 0.041847 0.94202 0.02696 0.1872
##Residuals 34   1.51039 0.044423         0.97304       
##Total     35   1.55223                  1.00000   

pshitreatpre10
##          Df SumsOfSqs  MeanSqs F.Model      R2   Pr(>F)   
##day        1    0.1658 0.165800  2.9582 0.08004 0.009998 **
##Residuals 34    1.9056 0.056048         0.91996            
##Total     35    2.0714                  1.00000    

pshitreatpre14
##         Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
##day        1   0.09563 0.095629   2.085 0.06302 0.0108 *
##Residuals 31   1.42183 0.045865         0.93698         
##Total     32   1.51746                  1.00000   

pshitreatpre28
##          Df SumsOfSqs  MeanSqs F.Model      R2  Pr(>F)  
##day        1   0.12542 0.125423  2.7528 0.11589 0.09375 .
##Residuals 21   0.95681 0.045562         0.88411          
##Total     22   1.08223                  1.00000   

pshitreatpre41
##          Df SumsOfSqs  MeanSqs F.Model      R2  Pr(>F)  
##day        1   0.10415 0.104147   2.315 0.08797 0.01172 *
##Residuals 24   1.07973 0.044989         0.91203          
##Total     25   1.18388                  1.00000     

pshitreat110
##          Df SumsOfSqs  MeanSqs F.Model      R2    Pr(>F)    
##day        1   0.20951 0.209515  3.8051 0.10065 0.0003999 ***
##Residuals 34   1.87208 0.055061         0.89935              
##Total     35   2.08160                  1.00000   

pshitreat1014
##          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
##day        1   0.14397 0.143969  2.5024 0.07469 0.0022 **
##Residuals 31   1.78353 0.057533         0.92531          
##Total     32   1.92750                  1.00000    

pshitreat1428
##          Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
##day        1   0.09558 0.095584  2.0612 0.10275 0.0625 .
##Residuals 18   0.83470 0.046372         0.89725         
##Total     19   0.93029                  1.00000   

pshitreat1441
##          Df SumsOfSqs  MeanSqs F.Model      R2  Pr(>F)  
##day        1   0.08669 0.086688   1.901 0.08301 0.01562 *
##Residuals 21   0.95762 0.045601         0.91699          
##Total     22   1.04431                  1.00000 

pshitreat2841
##         Df SumsOfSqs  MeanSqs F.Model      R2  Pr(>F)  
##day        1   0.08439 0.084391  1.8845 0.14626 0.03125 *
##Residuals 11   0.49260 0.044782         0.85374          
##Total     12   0.57699                  1.00000    


###### When worm-treated group changed, test for difference btwn pre28 and pre41 #####

#Make sure dfResil is set up with the correct levels. Should be inherited already, but just in case.
dfResil$day <- paste0("D", dfResil$day_num)
dfResil$day[dfResil$day == 'D-5'] <- 'pre'
dfResil$row <- paste(dfResil$mouse, dfResil$day, sep = "_") #Create the new row for merging with dfbetadiv.

#New dfResil already has beta diversity metrics bound to it, so can skip merging. 

#Ok, now subset the dfs with just the subgroups. 
dfResilhitreat <- dfResil[dfResil$microb=="HI" & dfResil$treatment=="worm",]
dfResilhitreat <- dfResilhitreat %>% tidyr::drop_na(ChangeASVs)
dfResilhictrl <- dfResil[dfResil$microb=="HI" & dfResil$treatment=="control",]
dfResilhictrl <- dfResilhictrl %>% tidyr::drop_na(ChangeASVs)
dfResillotreat <- dfResil[dfResil$microb=="LO" & dfResil$treatment=="worm",]
dfResillotreat <- dfResillotreat %>% tidyr::drop_na(ChangeASVs)
dfResilloctrl <- dfResil[dfResil$microb=="LO" & dfResil$treatment=="control",]
dfResilloctrl <- dfResilloctrl %>% tidyr::drop_na(ChangeASVs)
#
#Here, further subset the data. 
dfResillotreat2841 <- subset(dfResillotreat, day_num %in% c(28, 41))
dfResilhitreat2841 <- subset(dfResilhitreat, day_num %in% c(28, 41))
dfResilloctrl2841 <- subset(dfResilloctrl, day_num %in% c(28, 41))

#Bring all the 2841 datasets together for plotting. 
dfResil2841 <- rbind(dfResillotreat2841, dfResilloctrl2841, dfResilhitreat2841)
nrow(dfResil2841)
#
# No use making the hictrl version since there are no D28 HI ctrl.
View(dfResillotreat2841) #Check it out. 
#
#
#
#   ##### Now run the test to see whether there is a bigger difference in pre-28 or pre-41 alpha and beta diversity. #####
#
#      ##### Wilcoxon tests (did not do this for manuscript) #####
#LO worm group
wilcox.test(ChangeASVs ~ day_num, data=dfResillotreat2841)
##W = 192.5, p-value = 0.5599
wilcox.test(ChangeEvenness ~ day_num, data=dfResillotreat2841)
##W = 185, p-value = 0.7233
wilcox.test(Bray ~ day_num, data=dfResillotreat2841)
##W = 218, p-value = 0.179
wilcox.test(GUnifrac50 ~ day_num, data=dfResillotreat2841)
##W = 201, p-value = 0.4077


#HI worm group
wilcox.test(ChangeASVs ~ day_num, data=dfResilhitreat2841)
##W = 9.5, p-value = 0.1427
mod <- lme(GUnifrac50 ~ day, random=~1|mouse, 
           data = dfResilhitreat2841) 
summary(mod)
#For changeASVs
##                Value Std.Error DF   t-value p-value
##(Intercept) -42.86924 14.615153  7 -2.933205  0.0219
##dayD41       30.49424  9.325951  4  3.269826  0.0308
#For GUniFrac50
##                 Value  Std.Error DF   t-value p-value
##(Intercept) 0.27972406 0.02075590  7 13.476847  0.0000
##dayD41      0.00995135 0.02354976  4  0.422567  0.6943


wilcox.test(ChangeEvenness ~ day_num, data=dfResilhitreat2841)
##W = 26, p-value = 0.4351
wilcox.test(Bray ~ day_num, data=dfResilhitreat2841)
##W = 15, p-value = 0.5237
wilcox.test(GUnifrac50 ~ day_num, data=dfResilhitreat2841)
##W = 17, p-value = 0.7242


#LO sham group
wilcox.test(ChangeASVs ~ day_num, data=dfResilloctrl2841)
##W = 13, p-value = 0.4696
wilcox.test(ChangeEvenness ~ day_num, data=dfResilloctrl2841)
##W = 8, p-value = 0.132
wilcox.test(Bray ~ day_num, data=dfResilloctrl2841)
##W = 5, p-value = 0.04113
wilcox.test(GUnifrac50 ~ day_num, data=dfResilloctrl2841)
##W = 7, p-value = 0.09307



#
#

#      ##### Alternatively could run this as lmem (did this for manuscript) #####

#         ##### worm-treated LO #####
##dfResillotreat2841

#Copies 16S per NG gDNA
tmp16S <- dfResillotreat2841 %>% tidyr::drop_na(ChangeCopies16SperNG)
mod <- lme(ChangeCopies16SperNG ~ day, random=~1|mouse, data = tmp16S) 
summary(mod)
##               Value Std.Error DF   t-value p-value
##(Intercept) -2666979  376821.5 21 -7.077566       0
##dayD41       3078197  250059.0 13 12.309882       0

#ASV richness
mod <- lme(ChangeASVs ~ day, random=~1|mouse, data = dfResillotreat2841) 
summary(mod)
##               Value Std.Error DF   t-value p-value
##(Intercept) 37.49007  2.565653 22 14.612293  0.0000
##dayD41      -3.05529  2.499900 14 -1.222164  0.2418

mod <- lme(ChangeGUniFrac ~ day, random=~1|mouse, data = dfResillotreat2841) 
summary(mod)
##                 Value  Std.Error DF   t-value p-value
##(Intercept)  0.4251851 0.02257521 22 18.834159  0.0000
##dayD41      -0.0506774 0.01800933 14 -2.813953  0.0138


mod <- lme(ChangeEvenness ~ day, random=~1|mouse, data = dfResillotreat2841) 
summary(mod)
##                  Value  Std.Error DF    t-value p-value
##(Intercept)  0.01017553 0.02289691 22  0.4444061  0.6611
##dayD41      -0.03428656 0.02507271 14 -1.3674854  0.1930


mod <- lme(ChangeBray ~ day, random=~1|mouse, data = dfResillotreat2841) 
summary(mod)
##                 Value  Std.Error DF   t-value p-value
##(Intercept)  0.4513775 0.02773864 22 16.272517  0.0000
##dayD41      -0.0608215 0.02071164 14 -2.936583  0.0108

#         ##### worm-treated HI #####

tmp16S <- dfResilhitreat2841 %>% tidyr::drop_na(ChangeCopies16SperNG)
mod <- lme(ChangeCopies16SperNG ~ day, random=~1|mouse, data = tmp16S) 
summary(mod)
##               Value Std.Error DF    t-value p-value
##(Intercept)  304475.9  317537.9  6  0.9588647  0.3747
##dayD41      -296114.4  153333.6  1 -1.9311773  0.3042

##dfResilhitreat2841
mod <- lme(ChangeASVs ~ day, random=~1|mouse, data = dfResilhitreat2841) 
summary(mod)
##                Value Std.Error DF   t-value p-value
##(Intercept) -42.86924 14.615153  7 -2.933205  0.0219
##dayD41       30.49424  9.325951  4  3.269826  0.0308

mod <- lme(ChangeGUniFrac ~ day, random=~1|mouse, data = dfResilhitreat2841) 
summary(mod)
##                 Value  Std.Error DF   t-value p-value
##(Intercept) 0.27972406 0.02075590  7 13.476847  0.0000
##dayD41      0.00995135 0.02354976  4  0.422567  0.6943


mod <- lme(ChangeEvenness ~ day, random=~1|mouse, data = dfResilhitreat2841) 
summary(mod)
##                 Value  Std.Error DF    t-value p-value
##(Intercept)  0.01286894 0.01442265  7  0.8922729  0.4019
##dayD41      -0.01286385 0.01004742  4 -1.2803136  0.2696


mod <- lme(ChangeBray ~ day, random=~1|mouse, data = dfResilhitreat2841) 
summary(mod)
##                Value  Std.Error DF   t-value p-value
##(Intercept) 0.3250581 0.02319013  7 14.017092   0.000
##dayD41      0.0084554 0.02466783  4  0.342769   0.749


#         ##### sham-treated LO #####

##dfResilloctrl2841
tmp16S <- dfResilloctrl2841 %>% tidyr::drop_na(ChangeCopies16SperNG)
mod <- lme(ChangeCopies16SperNG ~ day, random=~1|mouse, data = tmp16S) 
summary(mod)
##               Value Std.Error DF   t-value p-value
##(Intercept) -1921418  297608.2  5 -6.456200  0.0013
##dayD41       1089427  297456.9  5  3.662471  0.0146

mod <- lme(ChangeASVs ~ day, random=~1|mouse, data = dfResilloctrl2841) 
summary(mod)
##                Value Std.Error DF   t-value p-value
##(Intercept) 12.500000  5.027137  5 2.486505  0.0554
##dayD41       5.666667  2.654145  5 2.135026  0.0859

mod <- lme(ChangeGUniFrac ~ day, random=~1|mouse, data = dfResilloctrl2841) 
summary(mod)
##                 Value  Std.Error DF   t-value p-value
##(Intercept) 0.27508947 0.02365699  5 11.628255  0.0001
##dayD41      0.06799771 0.02659196  5  2.557078  0.0508



mod <- lme(ChangeEvenness ~ day, random=~1|mouse, data = dfResilloctrl2841) 
summary(mod)
##                 Value  Std.Error DF    t-value p-value
##(Intercept) 0.01928692 0.02298736  5 0.8390231  0.4397
##dayD41      0.06813635 0.02251034  5 3.0268918  0.0292

mod <- lme(ChangeBray ~ day, random=~1|mouse, data = dfResilloctrl2841) 
summary(mod)
##                 Value  Std.Error DF  t-value p-value
##(Intercept) 0.28783333 0.02647495  5 10.87191  0.0001
##dayD41      0.08755556 0.03744123  5  2.33848  0.0665

#

########### Plotting for MS ############

#Create the MTD column. 
dfResil2841$MTD <- paste(dfResil2841$MT, dfResil2841$day, sep = "_")
dfResil2841$MTD <- factor(dfResil2841$MTD, levels=c("LO_worm_D28" ,"LO_worm_D41","LO_control_D28", "LO_control_D41", "HI_worm_D28","HI_worm_D41"))

#Boxplot change in richness pre-D28 and pre-D41
pResil2841ASVs <- ggplot(dfResil2841, aes(x=MTD, y=ChangeASVs)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfResil2841$ChangeASVs), max(dfResil2841$ChangeASVs)+0.1*max(dfResil2841$ChangeASVs)) +
  ylab("Change in ASV richness") +
  xlab("Group") +
  ggtitle("Change in ASV richness pre-28 and pre-41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  facet_wrap(~MT, scales="free_x")
pResil2841ASVs
#

#Boxplot change in qPCR 16S copies/ng DNA
dftmp <- dfResil2841 %>% tidyr::drop_na(ChangeCopies16SperNG)#Need to do this to use "min" and "max" functions in ylim
dftmp1 <- dfD41resist %>% tidyr::drop_na(ChangeCopies16SperNG)

pResil2841qPCR <- ggplot(dfResil2841, aes(x=MTD, y=(ChangeCopies16SperNG)/1000000)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
#  ylim(min(dfResil2841$ChangeCopies16SperNG), max(dfResil2841$ChangeCopies16SperNG)+0.1*max(dfResil2841$ChangeCopies16SperNG)) +
  ylim(min(dftmp$ChangeCopies16SperNG/1000000), (max(dftmp1$ChangeCopies16SperNG)/1000000) + (0.1*max(dftmp1$ChangeCopies16SperNG)/1000000)) +
  ylab("Change in 16S copies / ng DNA (in millions)") +
  xlab("Group") +
  ggtitle("Change in Copies16S per NG pre-28 and pre-41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  facet_wrap(~MT, scales="free_x")
pResil2841qPCR
#

#Boxplot change in richness pre-D28 and pre-D41
pResil2841GUnif <- ggplot(dfResil2841, aes(x=MTD, y=ChangeGUniFrac)) +
  geom_boxplot(aes(colour=MT), outlier.shape=NA) +
  geom_beeswarm(aes(color=MT, fill=MT, shape=MT), cex=2, size=3, alpha=0.9) +
  scale_shape_manual(values=c(16, 21, 15, 22)) +
  scale_colour_manual(values=c("#C66EF5","#C66EF5", "#FDD964" )) +
  scale_fill_manual(values=c("#C66EF5", "#FFFFFF", "#FDD964")) +
  #facet_wrap(~microb.x) +
  theme_bw() +
  ylim(min(dfD41resist$GU50), max(dfResil2841$ChangeGUniFrac)+0.1*max(dfResil2841$ChangeGUniFrac)) +
  ylab("Change in generalized UniFrac Distance") +
  xlab("Group") +
  ggtitle("Change in GUniFrac pre-28 and pre-41") +
  theme(plot.title=element_text(size=12, face="bold"), 
        legend.position= "none", 
        axis.title=element_text(size=14, face="bold"), 
        axis.text=element_text(size=13)) +
  facet_wrap(~MT, scales="free_x")
pResil2841GUnif
#

#      ##### Save the plots. #####
ggsave("~/Documents/1U4U_16S/GG2_reclassification/AlphaDiv/Resilience/pre28pre41_ASVchange_boxplots.pdf", pResil2841ASVs, units="in", height=5, width = 4)
ggsave("~/Documents/1U4U_qPCR/Resilience/pre28pre41_ChangeCopies16SperNG_boxplots.pdf", pResil2841qPCR, units="in", height=5, width = 4)
ggsave("~/Documents/1U4U_16S/GG2_reclassification/BetaDiv/Resilience/pre28pre41_GUniFracChange_boxplots.pdf", pResil2841GUnif, units="in", height=5, width = 4)
#
