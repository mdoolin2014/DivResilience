#Modifying my qiime script for performing it locally on my computer.

#Running with qiime2 on my own computer on 17Oct2021
#M.Doolin

#Activate qiime environment. Making this the same version used on the CHPC
# so that I know I'm getting the commands all right.
conda activate conda
source qiime2-2021.4
#or
conda activate qiime2-2021.4


#Which qiime environments can I have installed?
conda info --envs

#If you need to make a new classifier with a different taxonomy database, just
# get that over with before you start running your script. See my script in
# 16S_classifiers > MakeaQIIMEClassifier.sh


#######
#Ok now to get back to the processing.
cd ~/Documents/1U4U_16S/

mkdir GTDB_rerun1

#Define your common paths for the script
WRKDIR=/Users/mdoolin/Documents/1U4U_16S/GTDB_rerun
CLASSIFIER=/Users/mdoolin/Documents/16S_classifiers/GTDB_bac120_classifier_V3V4_14July22.qza
MANIFEST=/Users/mdoolin/Documents/1U4U_16S/Manifest_1U4U_Rerun_local.csv
METADATbyGNomexID=/Users/mdoolin/Documents/1U4U_16S/Rerun_metadata/MapByGNomexID_1U4U_w18684.tsv
METADATbySampleID=/Users/mdoolin/Documents/1U4U_16S/Rerun_metadata/MapBySampleID_1U4U_w18684.tsv


cd ${WRKDIR}

##### 1. IMPORT YOUR DATA
#Beware of which file type you're using. May give you an error for either
# tab-delimited or with extra spaces.
qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path  ${MANIFEST} \
 --output-path seqs_import.qza \
 --input-format PairedEndFastqManifestPhred33

# Visualize the reads
qiime demux summarize \
 --i-data seqs_import.qza \
 --o-visualization seqs_import.qzv



##### 2. REMOVE PRIMERS
#Run cutadapt on demuxed reads
qiime cutadapt trim-paired \
   --i-demultiplexed-sequences seqs_import.qza \
   --o-trimmed-sequences seqs_trim.qza \
   --p-front-f TGCCTACGGGNBGCASCAG \
   --p-front-r GCGACTACNVGGGTATCTAATCC \
   --p-cores 2 \
   --p-error-rate 0



##### 3. JOIN PAIRED END SEQUENCES
#Joining pairs. we're allowing maxdiffs to be pretty high bc of long seqs.
qiime vsearch join-pairs \
  --i-demultiplexed-seqs seqs_trim.qza \
  --o-joined-sequences seqs_trim_join.qza \
  --p-minmergelen 289 \
  --p-minovlen 20 \
  --p-maxdiffs 10 \
  --p-allowmergestagger \
  --verbose

#Visualize again to see what the reads look like
qiime demux summarize \
  --i-data seqs_trim_join.qza \
  --o-visualization seqs_trim_join.qzv
#Ok cool, so samples and merger looks great. No more reverse reads since there
# are only merged reads now.



##### 4. QUALITY FILTER AFTER JOINING.
#We will use deblur to denoise reads instead of DADA2. For deblur, you join the
# seqs first, and then denoise. In DADA2, it's the opposite.
qiime quality-filter q-score \
 --i-demux seqs_trim_join.qza \
 --p-min-quality 10 \
 --o-filtered-sequences seqs_trim_join_filt.qza \
 --o-filter-stats seqs_trim_join_filt_stats.qza

#view your q-score output
qiime metadata tabulate \
  --m-input-file seqs_trim_join_filt_stats.qza \
  --o-visualization seqs_trim_join_filt_stats-visualize.qzv

qiime tools view seqs_trim_join_filt_stats-visualize.qzv



##### 5. DENOISE WITH DEBLUR.
#If you're doing just V4, think about trimming to 250 bp, but for V3-V4, trim to
# 400-ish. Zac seems to trim to 392 for his things. Do this.
qiime deblur denoise-16S \
  --i-demultiplexed-seqs seqs_trim_join_filt.qza \
  --p-trim-length 392 \
  --p-jobs-to-start 4 \
  --o-table table.qza \
  --o-representative-sequences repseq.qza \
  --o-stats table_stats.qza

#Changing table.qza and repseq.qza to .qzv see the sequences of the strains.
#Great, that worked.
qiime feature-table summarize \
 --i-table table.qza \
 --o-visualization table.qzv

#Create a cool sequence table of the unique ASVs.
qiime feature-table tabulate-seqs \
--i-data repseq.qza \
--o-visualization repseq.qzv



##### 6. FILTER OUT CHIMERAS
#Filter out chimeras from your first table. Using methods to keep borderline chimeras.
qiime vsearch uchime-denovo \
--i-table table.qza \
--i-sequences repseq.qza \
--output-dir uchime_out

qiime feature-table filter-features \
--i-table table.qza \
--m-metadata-file uchime_out/chimeras.qza \
--p-exclude-ids \
--o-filtered-table table_nochim.qza

qiime feature-table filter-seqs \
--i-data repseq.qza \
--m-metadata-file uchime_out/chimeras.qza \
--p-exclude-ids \
--o-filtered-data repseq_nochim.qza

qiime feature-table summarize \
--i-table table_nochim.qza \
--o-visualization table_nochim.qzv

qiime feature-table tabulate-seqs \
--i-data repseq_nochim.qza \
--o-visualization repseq_nochim.qzv



##### 7. BUILD A PHYLOGENY
#Now, build phylogeny. This is important for doing any UniFRAC calculations
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences repseq_nochim.qza \
--o-alignment aligned_repseq.qza \
--o-masked-alignment masked_aligned_repseq.qza \
--o-tree tree_unroot.qza \
--o-rooted-tree tree_root.qza
#Here we made both a rooted and unrooted tree for these taxa.



##### 8. CLASSIFY YOUR READS BY A CLASSIFIER YOU'VE TRAINED.
#Ok now we're actually going to classify by the CLASSIFIER
qiime feature-classifier classify-sklearn \
--i-classifier ${CLASSIFIER} \
--i-reads repseq_nochim.qza \
--o-classification taxonomy.qza \
--p-n-jobs 2

#Make a visualizer after classification
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv


##### 9. FURTHER FILTERING OF YOUR TABLE
#Filter out any features (aka ASVs) that are singletons or doubletons.
qiime feature-table filter-features \
  --i-table table_nochim.qza \
  --p-min-samples 3 \
  --p-min-frequency 10 \
  --o-filtered-table doubleton-lowfreq-filt-table.qza
#Added 2 commands here, so they much have at least 3 total reads and be seen at
# least a total of 10 times across all samples.

#Filter out mitochondria and chloroplast
qiime taxa filter-table \
  --i-table doubleton-lowfreq-filt-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table final-table-byGNomexID.qza



##### 10. RENAME SAMPLES BASED ON METADATA COLUMN.
#Rename samples based on sample id instead of run ID for future labeling/figures.
#Will then have to make sure that calling samples in the future is done by new
# sample name.
qiime feature-table rename-ids \
  --i-table final-table-byGNomexID.qza \
  --m-metadata-file ${METADATbyGNomexID} \
  --m-metadata-column sample \
  --o-renamed-table final-table.qza
  #Note: If adding metadata in qiime2, it doesn't like using "yes" or "no", so make sure
  # that you don't use only those words as a binary category.



qiime metadata tabulate \
  --m-input-file final-table.qza \
  --o-visualization final-table.qzv


rm doubleton-lowfreq-filt-table.qza
rm table.qza


#Move to R, import into R with qiime2R package in R.


##### Visualizations and diversity investigation. #####

#Make a stacked bar.
qiime taxa barplot \
  --i-table final-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file ${METADATbySampleID} \
  --o-visualization table_taxbarplot.qzv


#To get core diversity metrics from this feature table
#This will output alpha and beta diversity metrics all into a folder for you.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny tree_root.qza \
  --i-table final-table.qza \
  --p-sampling-depth 3000 \
  --m-metadata-file ${METADATbySampleID} \
  --output-dir core-metrics-results-3000reads
#Choose sampling depth based on what you want to keep in.

#Pairwise beta diversity comparisons.
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-3000reads/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ${METADATbySampleID} \
  --m-metadata-column day \
  --o-visualization core-metrics-results-3000reads/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-3000reads/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file ${METADATbySampleID} \
  --m-metadata-column day \
  --o-visualization core-metrics-results-3000reads/weighted-unifrac-body-site-significance.qzv \
  --p-pairwise



##### Rarefying and subsetting dataset:

#Make sure you're in the right working directory.
cd ~/Documents/Anthelm_16S/qiime2_GTDB

#Rarefy your artefact.
qiime feature-table rarefy \
  --i-table final-table.qza \
  --p-sampling-depth 9174 \
  --o-rarefied-table psrare-table.qza

#If you want to look at it.
qiime metadata tabulate \
   --m-input-file qiime_outputs/psrare-table.qza \
   --o-visualization qiime_outputs/psrare-table.qzv

#Now subset rarefied dataset based on metadata for the samples that you want to keep.
# Just want to mimic psmtz phyloseq object, so making metadata with those samples.
# Remember must be a .tsv file.
qiime feature-table filter-samples \
  --i-table qiime_outputs/psrare-table.qza \
  --m-metadata-file ../metadata/psrareD7_metadat.tsv \
  --o-filtered-table qiime_outputs/psrareD7-table.qza

cd /Users/mdoolin/Documents/1U4U_16S/GTDB_rerun

#OR Subset your feature table based on a category in your whole.
qiime feature-table filter-samples \
   --i-table qiime_objects/psimptreat-table.qza \
   --m-metadata-file qiime_metadata/psexp_metadat.tsv \
   --p-where "[microb]='HI'" \
   --o-filtered-table qiime_objects/psimphitreat-table.qza

qiime feature-table filter-samples \
  --i-table qiime_objects/psimptreat-table.qza \
  --m-metadata-file qiime_metadata/psexp_metadat.tsv \
  --p-where "[microb]='LO'" \
  --o-filtered-table qiime_objects/psimplotreat-table.qza

#Control
qiime feature-table filter-samples \
  --i-table qiime_objects/psimpctrl-table.qza \
  --m-metadata-file qiime_metadata/psexp_metadat.tsv \
  --p-where "[microb]='HI'" \
  --o-filtered-table qiime_objects/psimphictrl-table.qza

qiime feature-table filter-samples \
  --i-table qiime_objects/psimpctrl-table.qza \
  --m-metadata-file qiime_metadata/psexp_metadat.tsv \
  --p-where "[microb]='LO'" \
  --o-filtered-table qiime_objects/psimploctrl-table.qza



#Check on your qiime object.
qiime metadata tabulate \
   --m-input-file qiime_objects/psimploctrl-table.qza \
   --o-visualization qiime_objects/psimploctrl-table.qzv

qiime feature-table summarize \
--i-table qiime_objects/psimp3000lotreat-table.qza \
--o-visualization qiime_objects/psimp3000lotreat-table.qzv










#
