#ANCOMBC in qiime2

#Running in qiime2-2022.11 using composition plugin

#set wd
cd /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification
METADAT=/Users/mdoolin/Documents/1U4U_16S/Rerun_metadata/psexp_metadat.txt

#Can run with one factor or multiple:
qiime composition ancombc --help
#if you want to see all options. 

#Example tables for this analysis. 
pspre1hitreat-nr
pspre1lotreat-nr
pspre1hictrl-nr
pspre1loctrl-nr

######### Run it for 1 file at a time. #########

#Multiple factors
qiime composition ancombc \
    --i-table qiime_objects/pairwise_tables/ps110lotreat-nr.qza \
    --m-metadata-file ${METADAT} \
    --p-formula day \
    --p-lib-cut 3000 \
    --p-p-adj-method 'holm' \
    --p-neg-lb \
    --o-differentials qiime_diffabund/ancombc_outputs/pslotreat110-nr-ancombc.qza
#The --p-neg-lb means that it includes structural 0's. If you instead had --p-no-neg-lb, it would drop structural 0's. But I want to include them since some important stuff
#    --p-neg-lb \
#

#Change the output to a spreadsheet type if you want to view through qiime viewer.
qiime composition tabulate \
 --i-data qiime_diffabund/ancombc_outputs/pslotreat110-nr-ancombc.qza \
 --o-visualization qiime_diffabund/ancombc_outputs/pslotreat110-nr-ancombc.qzv

#qiime tools view qiime_diffabund/ancombc_outputs/pslotreat110-nr-ancombc.qzv

#I'd like to see the spreadsheet of the results. I can export the files from this now. 
qiime tools export \
--input-path qiime_diffabund/ancombc_outputs/pslotreat110-nr-ancombc.qza \
--output-path qiime_diffabund/ancombc_outputs/pslotreat110-nr-ancombc

#
#


######### Write it in a loop ##########


#Try a loop to process a bunch of files without doing it manually. We've got a lot of pairwise tables.... 
cd /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification
METADAT=/Users/mdoolin/Documents/1U4U_16S/Rerun_metadata/psexp_metadat.txt

#Prep the output folder. 
mkdir qiime_diffabund/ancombc_outputs
#Run the ANCOM-BC. 
for table in qiime_objects/pairwise_tables/*nr.qza; do output="${table%}-ancombc.qza"
qiime composition ancombc \
  --i-table "$table" \
  --m-metadata-file ${METADAT} \
    --p-formula day \
    --p-lib-cut 3000 \
    --p-p-adj-method 'holm' \
    --p-neg-lb \
    --o-differentials "qiime_diffabund/ancombc_outputs/$output"
done

#Prep the output folder. 
mkdir qiime_diffabund/ancombc_outputs/exported_data
#Then write out the folders with the outcomes. 
for file in qiime_objects/pairwise_tables/*-ancombc.qza; do filename=$(basename "$file" .qza)
qiime tools export \
--input-path "$file" \
--output-path "qiime_diffabund/ancombc_outputs/exported_data/$filename"
done

#Now go to R to import files, concatenate, and merge with DESeq2 outputs. 
#
#
