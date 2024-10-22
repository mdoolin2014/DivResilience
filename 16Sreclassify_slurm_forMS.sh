#!/bin/bash
#SBATCH --account=penguin-np
#SBATCH --partition=penguin-shared-np
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120000
#SBATCH --time=8:00:00
#SBATCH -J 16S_reclassify
#SBATCH -o /uufs/chpc.utah.edu/common/home/dearing-group1/1U4U/qiime2/code/16S_reclassify.outerror

##### Submitted 17:55 on 4/28/23 job id 7685071 ####

#Reclassifying the taxonomy of all samples based on Greengenes2 instead of GTDB. Run on University of Utah CHPC.
#This script also inserts fragmented backbone for GG2 database that creates a better final phylogeny of all taxa. 

# Load qiime2. Source it in your miniconda3 env
module use ~/MyModules
module load miniconda3/latest
source activate qiime2-2022.11

# Set the scratch directory
SCRATCH=/scratch/general/vast/u1212572/Mar2023Run
WRKDIR=/uufs/chpc.utah.edu/common/home/dearing-group1/1U4U/reclassify_Apr2023
MANIFEST=${WRKDIR}/metadata/Manifest_1U4U_Rerun_slurm.csv
CLASSIFIERGG2=${WRKDIR}/2022.10.backbone.full-length_Prok34_classifier.qza
METADATbyGnomexID=/uufs/chpc.utah.edu/common/home/dearing-group1/1U4U/qiime2/metadata/MapByGNomexID_1U4U_w18684.tsv
METADATbySampleID=/uufs/chpc.utah.edu/common/home/dearing-group1/1U4U/qiime2/metadata/MapBySampleID_1U4U_w18684.tsv
RefTree=/uufs/chpc.utah.edu/common/home/round-group4/reference_seq_dbs/qiime2/ref_phylogenty_db/sepp-refs-gg-13-8.qza

cd ${WRKDIR}



##### 8. CLASSIFY YOUR READS BY A CLASSIFIER YOU'VE TRAINED.
#This is the classifier that Zac trained on greengenes2 database. 
#qiime feature-classifier classify-sklearn \
#  --i-classifier ${CLASSIFIERGG2} \
#  --i-reads repseq_nochim.qza \
#  --o-classification taxonomy_gg2.qza \
#  --p-n-jobs 8 \
#  --verbose

#Make a visualizer after classification
#qiime metadata tabulate \
#  --m-input-file taxonomy_gg2.qza \
#  --o-visualization taxonomy_gg2.qzv

# Taxa bar plots (no metadata map)
qiime taxa barplot \
 --i-table table_nochim.qza \
 --i-taxonomy taxonomy_gg2.qza \
 --m-metadata-file ${METADATbyGnomexID}
 --o-visualization table_tax_barplots_gg2.qzv


# Insertion fragment tree
qiime fragment-insertion sepp \
  --i-representative-sequences repseq_nochim.qza \
  --i-reference-database ${RefTree} \
  --o-tree tree_insertion_root.qza \
  --p-threads 8 \
  --o-placements insertion-placements.qza

qiime fragment-insertion filter-features \
  --i-table table_nochim.qza \
  --i-tree tree_insertion_root.qza \
  --o-filtered-table table_ff.qza \
  --verbose \
  --o-removed-table table_uninserted_frags.qza

qiime feature-table summarize \
  --i-table table_uninserted_frags.qza \
  --o-visualization table_uninserted_frags.qzv
#If these aren't any different, then table_nochim is the same as table_ff. 



##### 9. FURTHER FILTERING OF YOUR TABLE
#Filter out any features (aka ASVs) that are singletons or doubletons.
qiime feature-table filter-features \
    --i-table table_ff.qza \
    --p-min-samples 3 \
    --p-min-frequency 10 \
    --o-filtered-table doubleton-lowfreq-filt-table.qza
#Added 2 commands here, so they much have at least 3 total reads and be seen at
# least a total of 10 times across all samples.

#Filter out mitochondria and chloroplast reads
qiime taxa filter-table \
    --i-table doubleton-lowfreq-filt-table.qza \
    --i-taxonomy taxonomy_gg2.qza \
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



##### Visualizations and diversity investigation. #####

#Make a stacked bar.
qiime taxa barplot \
  --i-table final-table.qza \
  --i-taxonomy taxonomy_gg2.qza \
  --m-metadata-file ${METADATbySampleID} \
  --o-visualization final-table_taxbarplot.qzv


#To get core diversity metrics from this feature table
#This will output alpha and beta diversity metrics all into a folder for you.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny tree_insertion_root.qza \
  --i-table final-table.qza \
  --p-sampling-depth 3000 \
  --m-metadata-file ${METADATBySampleID} \
  --output-dir core-metrics-results
#Choose sampling depth based on what you want to keep in.


#Pairwise beta diversity comparisons.
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file ${METADATBySampleID} \
  --m-metadata-column day \
  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column day \
  --o-visualization core-metrics-results/weighted-unifrac-body-site-significance.qzv \
  --p-pairwise










##
##