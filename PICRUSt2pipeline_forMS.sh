#PICRUSt2 whole pipeline, not in QIIME2
#Created 29Mar22, from https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.4.2)
#Similar tutorial here: https://github.com/LangilleLab/microbiome_helper/wiki/CBW-2021-PICRUSt2-Tutorial
#Outside of qiime, just install with bioconda.
#Most recently run September 2024 for diversity/resilience project#

#need to add Python to PATH before starting because somehow my 
/Users/mdoolin/Library/Python/3.9/bin

echo export PATH="/Users/mdoolin/Library/Python/3.9/bin:$PATH" >> ~/.zshrc
#or can add it in manually by opening and editing the .zshrc file directly. 
nano .zshrc
#save and exit.


echo export PATH="<PATH_TO_PYTHON>:$PATH" >> ~/.profile

#Install picrust with mamba, from source. 
wget https://github.com/picrust/picrust2/archive/v2.5.3.tar.gz
tar xvzf  v2.5.3.tar.gz
cd picrust2-2.5.3/

mamba env create -f picrust2-env.yaml
conda activate picrust2
pip install --editable .
#Test install. Seems to search the whole computer, not sure if this is necessary...
pytest 



#First, you need to make your qiime outputs into the right file types.
#They're using a tutorial dataset, but I'm going to use my 1U4U outputs from qiime.
# For this, need metadata.tsv, seqs.fna, table.biom. Will just have to convert repseqs
# to fasta (.fna)
cd /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification/qiime_objects/

conda activate qiime2-2022.11
mkdir ./repseq-exports


qiime tools export \
  --input-path repseq_nochim.qza \
  --output-path repseq_nochim-exports
  #this exports dna-sequences.fasta

conda deactivate

#Activate picrust2 environment.
conda activate picrust2

cd /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification/biom_inputs

#Check out your input files to make sure they're set up correctly.
#Biom file needs to have ASVs as rows and samples as columns (for the full pipeline v2.5.3)
biom head -i psexp3000.biom
biom summarize-table -i pspre41lotreat-nr.biom
less ../qiime_objects/repseq_nochim-exports/dna-sequences.fasta

mkdir ../picrust2_out_pipeline
cd ../picrust2_out_pipeline



###If you want, you can just run the whole pipeline in one command to get default outputs:
picrust2_pipeline.py -s study_seqs.fna -i study_seqs.biom -o picrust2_out_pipeline -p 1



### Or you can do it step-by-step. 

#ok, now place your samples into a reference tree.
place_seqs.py -s ../qiime_objects/repseq_nochim-exports/dna-sequences.fasta \
-o out.tre \
-p 2 \
--intermediate intermediate/place_seqs
#Note that for repseq_nochim from 1U4U dataset, only 1 sequence did not align well to the reference. 

#Note that I had some issues with dependencies. uninstalled r-jsonlite and figured it out after several hours....


#Run hsp.py for hidden state prediction, i.e., to predict the missing genomes
# of each ASV, aka predict the number of gene families.
#Predict how many 16S copies there are associated with the predicted genome of each sequence
hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 2 -n
#And predict kegg orthologs.
hsp.py -i KO -t out.tre -o KO_predicted.tsv.gz -p 2
#And predict enzyme classifications (EC)
hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 2

#Check the outputs to see what they look like, did they come out ok?
zless -S marker_predicted_and_nsti.tsv.gz
zless -S KO_predicted.tsv.gz


mkdir picrust2_out_pipeline/psexp3000_metagenome_out


#These are actually usable outputs if you want to just look at counts.
# But if you want to do metagenome inferences, need to do more steps.

#Make sure the biom file is set up correctly and not transposed. Need samples as rows
# and ASV counts as columns.

metagenome_pipeline.py \
-i /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification/biom_inputs/psexp3000.biom \
-m marker_predicted_and_nsti.tsv.gz \
--strat_out \
-o /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification/picrust2_out_pipeline/psexp3000_metagenome_out \
-f KO_predicted.tsv.gz  
# --strat-out means a stratified output will be output (as contrib), in addition
# to the unstrat output.
#psexp3000 all ASVs were below the max NSTI cutoff of 2.0, so all were retained for downstream analyses. 


#Take a look at an output to make sure everything worked.
zless -S psexp3000_metagenome_out/pred_metagenome_unstrat.tsv.gz

#This is to create a legacy format table from PICRUSt1 so you can run BURRITO and
# MIMOSA on it.
convert_table.py psexp3000_metagenome_out/pred_metagenome_contrib.tsv.gz \
                 -c contrib_to_legacy \
                 -o pred_metagenome_contrib.legacy.tsv.gz

#Final step is to infer the pathway abundance from the metagenome predictions.
#Note that this is only for EC pathways, not KO. 
pathway_pipeline.py -i psexp3000_metagenome_out/pred_metagenome_contrib.tsv.gz \
                 -o pathways_out -p 1

#Combined outputs for the different picrust outputs into one folder – picrust_out.
#Making sure I'm now in that output folder as my working directory.
cd /Users/mdoolin/Documents/1U4U_16S/GG2_reclassification/picrust2_out_pipeline/psexp3000_metagenome_out

#create a path to the mapfile of your choice, if not using default
#mapfile = /Users/mdoolin/picrust2-2.5.3/picrust2/default_files/description_mapfiles

#Try -m KEGG_pathways_info.tsv.gz or call it locally after downloading with --custom_map_table from a directory of your choosing. This give the high-level ko identities that are more easily summarized rather than the specific pathway names. Could probably also add the regular pathway names, too, if you wanted both.  
#Ok and then just add descriptions for each functional ID in abundance tables.
add_descriptions.py -i pred_metagenome_unstrat.tsv.gz -m KO  \
                    -o pred_metagenome_unstrat_descrip.tsv.gz


add_descriptions.py -i path_abun_unstrat.tsv -m METACYC \
                    -o path_abun_unstrat_descrip.tsv




#Ok now you can transition these outputs to a phyloseq object in R. 
