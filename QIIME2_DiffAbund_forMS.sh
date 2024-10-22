#Differential abundance testing in qiime2 
#Running gneiss and ANCOM, although I think gneiss is defunct or soon to be so.
#See Moving Pictures tutorial. 

#Note if you want to export the data from your qiime object, use this:

cd /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification


#Note that I have created the script "QIIME_FeatureTableSubsetting.sh" for code to load in feature tables from BIOM files in R and to subset and view feature tables. 


#Run gneiss clustering tutorial to get an idea of the relationships of samples.
https://docs.qiime2.org/2022.8/tutorials/gneiss/
#You can either cluster with the correlation clustering. This is just generally
# clustering them, not by any variable in particular:
qiime gneiss correlation-clustering \
  --i-table GG2_reclassification/qiime_objects/psexp3000hitreat.qza \
  --o-clustering GG2_reclassification/qiime_diversity/psexp3000hitreat-cluster.qza

#Or with the gradient-clustering, to group samples based on some numeric category.
# Think something like age or day. Note from tutorial
#An important consideration for downstream analyses is the problem of overfitting.
# When using gradient-clustering, you are creating a tree to best highlight
# compositional differences along the metadata category of your choice, and it is
# possible to get false positives:
qiime gneiss gradient-clustering \
  --i-table GG2_reclassification/qiime_objects/psexp3000hitreat.qza \
  --m-gradient-file Rerun_metadata/psexp_metadat.tsv \
  --m-gradient-column day_num \
  --o-clustering  GG2_reclassification/qiime_diversity/psexp3000hitreat-gradient-cluster.qza \
  --p-no-weighted

#And then make the heatmap to visualize with whichever clustering mechanism you want.
#The metadata column needs to be categorical for this one.
qiime gneiss dendrogram-heatmap \
  --i-table GG2_reclassification/qiime_objects/psexp3000hitreat.qza \
  --i-tree GG2_reclassification/qiime_diversity/psexp3000hitreat-gradient-cluster.qza \
  --m-metadata-file Rerun_metadata/psexp_metadat.tsv \
  --m-metadata-column infection_status \
  --p-color-map seismic \
  --o-visualization GG2_reclassification/qiime_diversity/psexp3000hitreat-gradientcluster-heatmap.qzv

qiime tools view qiime_diversity/psexp3000hitreat-heatmap.qzv




#
#
#Run ANCOM from Moving Pictures tutorial.
https://docs.qiime2.org/2022.11/tutorials/moving-pictures/

#Look at the table you're using to make sure it's got the samples you want. 
qiime metadata tabulate \
   --m-input-file biom_inputs/pspre14hictrl-nr.qza \
   --o-visualization qiime_objects/pspre14hictrl-nr.qzv

#Still in ~/Documents/Anthelm_16S/qiime2_GTDB
cd ~/Documents/1U4U_16S/GG2_reclassification
#First, need to convert feature table into compositional version from the count table.
qiime composition add-pseudocount \
  --i-table qiime_objects/pspre14lotreat-nr.qza \
  --o-composition-table qiime_diffabund/pspre14lotreat-comp-table.qza

#Metadata column here also needs to be categorical.
qiime composition ancom \
  --i-table qiime_diffabund/pspre14loctrl-comp-table.qza \
  --m-metadata-file ../Rerun_metadata/psexp_metadat.tsv \
  --m-metadata-column day \
  --o-visualization qiime_diffabund/pspre14loctrl-ancom.qzv

qiime tools view diffabund/pspre41hitreat-ancom.qzv



#Or run ANCOM at a different taxonomic level.
cd ~/Documents/1U4U_16S/GG2_reclassification
#Level 2 is phylum, 5 is family, 6 is genus, 7 is species
qiime taxa collapse \
  --i-table qiime_objects/pspre41hictrl-nr.qza \
  --i-taxonomy qiime_objects/taxonomy_gg2.qza \
  --p-level 5 \
  --o-collapsed-table qiime_objects/pspre41hictrl-nr-family-table.qza

qiime composition add-pseudocount \
  --i-table qiime_objects/pspre41hictrl-nr-family-table.qza \
  --o-composition-table qiime_diffabund/pspre41hictrl-nr-family-comp-table.qza

qiime composition ancom \
  --i-table qiime_diffabund/pspre41lotreat-nr-family-comp-table.qza \
  --m-metadata-file ../Rerun_metadata/psexp_metadat.tsv \
  --m-metadata-column day \
  --o-visualization qiime_diffabund/pspre41lotreat-nr-family-ancom.qzv

qiime tools view qiime_diversity/diff_abund/pspre41hitreat-nr-family-ancom.qzv


#
#

#Coming back to this on 4 April 2024
#Or run ANCOMBC rather than original ANCOM. The BC stands for "Bias Correction"
#Running in qiime2-2022.11 using composition plugin

#set wd
cd ~/Documents/1U4U_16S

#Can run with one factor or multiple:
qiime composition ancombc --help
#if you want to see all options. 

#Multiple factors
qiime composition ancombc \
    --i-table GG2_reclassification/qiime_objects/psexp-control-table.qza \
    --m-metadata-file Rerun_metadata/psexp_metadat.tsv \
    --p-formula 'microb + infection_status' \
    --p-reference-levels microb::HI \
    --o-differentials GG2_reclassification/qiime_diffabund/psexp-control-ANCOMBC.qza


pspre41hitreat-nr
pspre41lotreat-nr
pspre41hictrl-nr
pspre41loctrl-nr


#Single factor, this is how I'm testing two timepoints within a microbiome/treatment group.
qiime composition ancombc \
    --i-table GG2_reclassification/qiime_objects/pspre41loctrl-nr.qza \
    --m-metadata-file Rerun_metadata/psexp_metadat.tsv \
    --p-formula day \
    --o-differentials GG2_reclassification/qiime_diffabund/pspre41loctrl-nr-ANCOMBC.qza


#Change the output to a spreadsheet type.
qiime composition tabulate \
 --i-data GG2_reclassification/qiime_diffabund/pspre41loctrl-nr-ANCOMBC.qza \
 --o-visualization GG2_reclassification/qiime_diffabund/pspre41loctrl-nr-ANCOMBCviz.qzv


#I'd like to see the spreadsheet of the results. I can export the files from this now. 
qiime tools export \
--input-path GG2_reclassification/qiime_diffabund/pspre41loctrl-nr-ANCOMBC.qza \
--output-path GG2_reclassification/qiime_diffabund/pspre41loctrl-nr-ANCOMBCtable

#
