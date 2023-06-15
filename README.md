# Evolutionary_history_CaSR
Evolutionary history of Calcium-sensing receptors sheds light into hyper/hypocalcemia-causing mutations
## Retrieve all class C GPCRs

Proteomes were downloaded from [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/)


To get all class C GPCRs 7tm_3 domain profile was used. Each species proteome was searched against `7tm_3.hmm`. Following python code takes each species proteome file, performed hmmseach and saves the result. For example, the human proteome with the ID `9606` was searched against `7tm_3.hmm`, and the result file was saved as *9606.out*.
 ```
 python do_hmmsearch.py
 
 ```
 To retrieve sequences that match the 7tm_3.hmm (above the default HMM threshold), Python code was used to parse the hmmsearch result files. Each output was saved as speciesID.fasta.
 
 ```
  python get_fasta_files_for_hmmscan.py 
 ```

Pfam-A profile was used to obtain domain architectures of class C GPCRs that were retrieved by using *get_fasta_files_for_hmmscan.py* code. Each *speciesID.fasta* file was scanned against *Pfam-A.hmm*

 ```
 python do_hmmer_scan.py
 ```

Since NCBI proteomes include isoforms, longest isoforms were selected as canonicals. For Human we selected canonical sequences given in Uniprot. 

 ```
 python take_single_isoform.py
 ```
 
We added taxid and geneid to fasta files. Based on hmmscan results, we parsed the 7TM domains of sequences to use later.

```
python get_geneid_taxid.py
python get_7tm_3_domain.py
 ```
For each species a blast database was built by using single isoforms.

 ```
makeblastdb -in speciesid_single_isoforms.fasta -parse_seqids -dbtype prot -out speciesid
 ```
 
Each sequence of each species was quired through BLASTP against each species blast database. 

 ```
python do_blastp.py queryID subjectID
 ```
From blast results bi-directional mutual best hits of each human class C GPCR subfamily were retrieved.

 ```
python retrieve_mutual_best_hits.py
 ```
 
Since a species was not blasted against its own blastdatabase, canonical human GPCR was added to its bi-directional mutual best hits. 7TM domains of bi-directional mutual best hits of each human class C GPCR subfamily were trimmed.
 ```
python get_7tm3_domain_of_mutual_hits.py
 ```
7TM domains of each subfamily were aligned with *mafft einsi* and ML trees were built by the following command:
 ```
raxmlHPC -s msa_file -w output_path -m PROTGAMMAAUTO -p 12345 -x 12345 -c 4 -#10 --no-seq-check -f a -n output_file 
 ```
Common ancestor of each node was added to the ML trees by using `internal_node_info.py` code.

Full length sequences of each subfamily were aligned and trimmed to build a HMM profile. 
 ```
python build_hmmer_profiles.py
 ```
Class C single isoforms were compiled into a single fasta file and hmmscan was perfomed against subfamily profile HMMs. HMMscan results were saved into a json file.
 ```
python profile_hmmscan.py
python make_hmmscan_json_file.py

 ```
 
Each sequence was assigned to a subfamily (1) The maximum score value of the hmmscan belonged to the given subfamily, (2)E-value is a measure of the significance of a match in a database search and the lower the E-value, the more significant the match is. The E-value of the sequence must be the lowest. (3)The sequence must belong to the most common highest taxonomic level of the given subfamily. 


```
python retrieve_subclass_sequences.py subfamily_tree_file subfamily_name
```

To remove paralogs from the subfamilies, all retrieved sequences were aligned with MAFFT `einsi`. ML tree was built with the following command:
```
raxml-ng --msa subfamily.fasta data-type AA --model JTT --prefix subfamily --seed 2 --thread 8
```

Duplication nodes were labeled with `label_duplication_nodes` code. From duplication  labeled trees, trees were pruned by using FigTree and Ete3. 

After paralog removal,  264 CaSR_like sequences and all GPRC6A, CaSR and TAS1R1-3 proteins were aligned and ML tree was built:

```
iqtree2 -s aln.fasta -m MFP -bb 1000 -nt AUTO
```

## Predict the Mutation Types in CaSR 
To predict the mutation type (LoF/GoF) in human CaSR,a gradient boosting- based machine learning algorithm, XGBoost was used.

For each dataset split the sklearn train test split model was used with stratify option to keep loss-of-function to gain-of- function ratio almost the same in the datasets.

Subfamily alignments and mutations were divided randomly as 80% training and the remaining 20% test data before creating feature matrices to prevent information leakage because conservation scores were used as features to train our model. 25% of the training data was randomly picked as the validation data five times for cross validation. 

To create 50 different csv files:

```
python create_ml_csv.py 
```

To train model, *run_xgb_varImp_feature_elimination.R* code was used by giving id. To train all splitted datasets, the code was called 50 times by giving ids 1 to 50.  

```
Rscript run_xgb_varImp_feature_elimination.R $id 
```

Results were saved as Rdata files: casr_xgb_replication_id_prauc.RData, casr_xgb_replication_id_prediction.RData,casr_xgb_replication_id_result.RData,casr_xgb_replication_id_state.RData,casr_xgb_replication_id_xgb.model 

To plot AUROC and AUPR for train and test sets, Rdata files were also saved as csv with Results_to_csv.R and organized the csv files by usign *casr_ml_parse_auc.py* for simplicity.

`casr_ml_results.R` code draws AUROC and AUPR graphs for train and test sets.

After model performance was reportes, 5 fold-cross validation was used to select features to train whole dataset.

*create_final_test_csv.py* creates the 5 validation and one whole dataset csv files. Final features were selected by running:
```
Rscript feature_selection_grid_search_cv.R 1
```

To make predictions for all substitutions at each position, *neutral* substitutions were eliminated and final test csv was created from a new UNIPROT alignments by usin *create_final_test_csv.py* code. 

Whole dataset trained and test predictions were made by *train_whole_data.R* code. Whole train and test predictions scores were saved as csv files.

To find feature importance SHAP values were calculated and graph was plotted by *shap_values.R*

Test and whole train csv files were also created without BLOSUM encoding to show amino acid names by `dummy_aa_csv.py` and `dummy_aa_test_csv.py` codes.

To draw heatmap, `R_heatmap_csv.py` code was used to prepare a csv file includes 'position', 'substitution', 'prediction','result'. *test_heatmap.R* code draws the heatmap and substituton count graphs. 



## SDP Scores

To calculate SDP (Site-wise Differential Profile) scores, the alignment of 264 CaSR-like sequences, including CaSR, GPRC6A, and TAS1R1-3, and their ML (Maximum Likelihood) tree were used.

First, positions in the alignment that were correspond to gaps in the human sequence were removed. This modified alignment file was used for ancestral sequence reconstruction. An example command for reconstructing ancestral sequences for CaSR is as follows:

```
iqtree2 -s casr_6a_tastes_nogapcasr.fasta  fasta -te casr_group_iqtree_rerooted.nwk -m JTT+R10 -asr -pre ancestral_no_gap_casr
```

Using the ancestral tree obtained from the previous step (.treefile file), each subfamily tree was pruned and saved as a separate newick file.

In a text file, the node name of each subfamily and the distances between the target subfamily node and other subfamily nodes were recorded. The node names were automatically assigned by IQTree.

The Test_SpecificPos.R code was used to calculate SDP scores based on the recorded information.


## Subfamily Specific Profile HMMs

To make subfamily specific profile HMMs, firstly from the alignment of the subfamily, default HMM model must be produced with `--pnone` option of hmmer. Secondly, identity scores must be calculated with `compute_identity_score.py` code. Then, `Test_HMM.R` can be run to produce subfamily specific profile HMMs.