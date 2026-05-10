#Processing mouse blood serum metabolite data from GCMS

##### Import normalized and scaled data #####

setwd("/Users/mdoolin/Documents/1U4U_metab/GCMS")

gcnormall <- read.csv("normalizedGCMSoutputs.csv")
gcnorm <- subset(gcnormall, !label %in% c("pb", "qc"))
rownames(gcnorm) <- gcnorm$samlab  #Make the sample IDs into the rownames.
gcscaledall <- read.csv("scaledGCMSoutputs.csv")
gcscaled <- gcscaledall[gcscaledall$label != "qc",]
rownames(gcscaled) <- gcscaled$samlab  #Make the sample IDs into the rownames.


##### Create a PCA plot of GCMS data ##### 

cols <- c('#A11FE5', '#440b63', '#440b63', "#DCAE1C", "#6e5509", "#6e5509", '#000000')

pp <- autoplot(scaled.pca, data=gcscaled, fill="label", color="label", shape='microb',
               size=3, x=1, y=2) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=fills) +
  stat_ellipse(aes(colour=label), type="norm", level=0.9) +
  scale_shape_manual(values=c(22,21)) +
  ggtitle("All Normalized, scaled PCA") +
  theme_bw()
#mess around with plotting params https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html 

pp


#ggsave('plots/AllScaledGCMS_PCA_correctcolors.pdf', pp, units="in", height=6, width=8)


##### Example of composition analyses #####

D41exp.norm <- subset(gcnorm, treatment=="parasite" & timepoint=="cleared")

D41exp.mat <- as.matrix(D41exp.norm[,-c(1:9)])  #create a matrix from the df. Numeric only
D41exp.rel <- as.data.frame(funrar::make_relative(D41exp.mat)) #make that relative abundances. 
# Is relative by row, which is correct. 
D41exp.metadat <- D41exp[,c(1:9)] 
D41exp.rel <- cbind(D41exp.metadat, D41exp.rel)

bray <- vegdist(D41exp.rel[,-c(1:9)], method="bray")

#Normality testing for PERMANOVA with betadisper
treatments <- D41exp.rel$label
mod <- betadisper(man, treatments)
mod
plot(mod)


#PERMANOVA
perm <- adonis2(bray ~ label, data = D41exp.metadat)
perm


