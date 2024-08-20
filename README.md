# 28S_marker_paper
 This is a repository of the code and data used in the analysis of the manuscript entitled: "Introducing a novel 28S rRNA marker for improved assessment of coral reef biodiversity".

## Primer_Alignment
This folder contains the aligned sponge sequences used to design the novel primer pair.

## Fastq Sequences
The raw, demultiplexed Fastq files were submitted to NCBI. This project's NCBI Sequence Read Archive (SRA) numbers are SRR30275729-SRR30275738; SRR30316297-SRR30316304; SRR30316420-SRR30316422 under Bioproject PRJNA1149044 and can be obtained here: https://www.ncbi.nlm.nih.gov/sra.

Data were sequenced in three MiSeq runs: 1) 28S ARMS A and associated mock communities, 2) 28S ARMS B and C and associated mock communities together with *Mock community 1*, and 3) COI ARMS A, B and C data. The Fastq file names, corresponding sample names and treatments are in the Martineau_Testing_primers_sequence_metadata.csv file in the Bioinformatics folder. 

Sequences were processed using the DADA2 pipeline. See the manuscript for details. The rarefied, clean abundance tables' outputs and the BLASTn outputs are in the Bioinformatics folder. 


## Bioinformatics
The Bioinformatics folder is divided into 3 subfolders, each containing all the data and code needed to create the taxonomy tables for the 3 different sequencing runs.  Refer to manuscript for additional information.

### 1. ARMSA28S

contains the post-DADA2 abundance tables, BLASTn outputs, and custom R script to assign taxonomy, filter, and merge BLASTn results for ARMS A and associated mock communities. The resulting taxonomy and abundance tables are located in the 28SCOI folder. 

- Cleaned post-DADA2 abundance table- Martineau_TestingPrimers_Clean_ARM17_abundance_071224.rarefied.csv
- Cleaned post-DADA2 abundance table with sequences - Martineau_TestingPrimers_Clean_ARM17_abundance_sequences_071224.rarefied.csv
- BLASTn output of sequences queried against local sponge database - blast_results_ARMS17.tab
- BLASTn output of sequences queried against GenBank nt database - blast_results_ARMS17_nt.tab
- R script filtering BLASTn outputs, assigning taxonomy, merging results and duplicated taxonomic entities - Martineau_TestingPrimers_bioinformatics_ARMSA28S_070224.R

### 2. COI
contains the post-DADA2 abundance tables, BLASTn outputs, and custom R script to assign taxonomy, filter, and merge BLASTn results for COI data. The resulting taxonomy and abundance tables are located in the 28SCOI folder. 

- Cleaned post-DADA2 abundance table- Testing_primers_abundance_ASV_COI_070224.rarefied.csv
- Cleaned post-DADA2 abundance table with sequences - Testing_primers_sequence_ASV_COI_070224.rarefied.csv
- BLASTn output of sequences queried against local sponge database - blast_results_COI_sponges_070224.tab
- BLASTn output of sequences queried against GenBank nt database - blastn_COI_outputs-070224.tab
- R script filtering BLASTn outputs, assigning taxonomy, merging results and duplicated taxonomic entities - Martineau_TestingPrimers_bioinformatics_COI.R

### 3. 28SCOI
contains the post-DADA2 abundance tables, BLASTn outputs, and custom R script to assign taxonomy, filter, and merge BLASTn results for *Mock 	community 1*, 28S ARMS B and C, and associated mock communities. This R script also merges abundance and taxonomy tables from the 3 sequencing runs into one clean and final output with no duplicated taxonomic entities. Those final tables are located in the Clean_tables folder. 

- Cleaned post-DADA2 abundance table- Martineau_Testing_primers_abundance_ASV_070224.rarefied.csv
- Cleaned post-DADA2 abundance table with sequences - Martineau_Testing_primers_abundance_sequences_070224.rarefied.csv
- BLASTn output of sequences queried against local sponge database - blast_results_spongesdb.tab
- BLASTn output of sequences queried against GenBank nt database - blast_results_ntdb.tab
- Clean abundance and taxonomy tables for COI sequencing run - sponge_ntCOI_taxonomy_table.csv and sponge_ntCOI_abundance_table.csv
- Clean abundance and taxonomy tables for ARMS A and associated mock communities - sponge_ntARMS17_taxonomy_table.csv and sponge_ntARMS17_abundance_table.csv
- Meta data information about sample names: Martineau_TestingPrimers_28SANDCO1_Metadata.csv
- R script filtering BLASTn outputs, assigning taxonomy, merging results of 3 sequencing runs, and deleting duplicated taxonomic entities - Martineau_TestingPrimers_bioinformatics_28SCOI_070224.R


## Clean_tables
This folder contains all data (the clean taxonomy, abundance tables, and metadata files) used for manuscript data analyses.  

- Taxonomy, abundance tables, and metadata files for Figure S1: sponge_ntARMS17_taxonomy_table.csv, sponge_ntARMS17_abundance_table.csv and Martineau_TestingPrimers_metadata_17.csv
Taxonomy and abundance tables were produced in the Martineau_TestingPrimers_bioinformatics_ARMSA28S_070224.R script. 

- Taxonomy, abundance tables, and metadata files for Figure 2: taxonomy_all_OTUs_070224.csv, abundance_all_OTUs_070224.csv, and Martineau_TestingPrimers_CO1AND28S_metadata_allOTUs.csv. 
Taxonomy and abundance tables were produced in the Martineau_TestingPrimers_main_figures_070224.R script. 

- Taxonomy, abundance tables, and metadata files for all other figures and analyses: combined28SCOI_taxonomy_table_071224.csv, combined28SCOI_abundance_table_071224.csv, Martineau_TestingPrimers_CO1AND28S_metadata.csv.
Taxonomy and abundance tables were produced in the Martineau_TestingPrimers_bioinformatics_28SCOI_070224.R script. 



## Data_Analyses
All  R codes used to obtain the results are presented in the Data_analyses folder. All Figures and accompanying Tables presented in the manuscript are labelled in the R files and are found here:

- Figure 2: Martineau_TestingPrimers_Figure2_070224.R
- Figure S1: Martineau_TestingPrimers_FigureS1_071224.R
- Figures 3-6, S2-S3, Tables S6-S13: Martineau_TestingPrimers_main_figures_070224.R
