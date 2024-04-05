
# IAV ancestral reconstruction at the swine-human interface

Pipeline, scripts, and data for the paper:

_Anker, K.M., Ciucani, M.M., Nissen, J.N., Anderson, T.K., Pedersen, A.G., and Trebbien, R._ . **Exploring genetic signatures of zoonotic influenza A virus at the swine-human interface with phylogenetic and ancestral sequence reconstruction** (in preparation).



## Analysis pipeline and scripts
The steps of the pipeline for data preparation and analysis is described in the [ancestor_project_notes](https://github.com/KMAnker/IAV_ancestral_reconstruction/blob/main/ancestor_project_notes.sh) file.
Here, it can also be seen where each of the scripts in the [scripts](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/scripts) folder is used.



## Data
Sequence data for this analysis was obtained from GISAID and NCBI Influenza Virus Resource (Genbank).

We gratefully acknowledge all data contributors, i.e., the authors and their originating laboratories responsible for obtaining the specimens, and their submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID initiative as well as the NCBI Influenza Virus Resource database, on which this research is based.


We cannot publicly upload the genetic sequences from GISAID, and therefore no sequence (fasta) files are uploaded to this repository. However, the [data](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/data) folder contains lists of accession numbers for all the sequences originally downloaded for analysis.



## Analysis and output folders
Folders containing input and output files for the different steps of the analysis:
- [treefiles](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/treefiles): IQtree output files and both input and output files from Treetime analysis (except fasta sequence files).
- [anclib_files](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/anclib_files): Output files from the anclib analysis.
- [selection_analysis](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/selection_analysis): Output json files from the HyPhy aBSREL and MEME analyses as well as text and excel files containing concatenated output from all segments.
- [xgboost](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/xgboost): one-hot encoded input files and output files from xgboost analyses (run with the [xgboost.ipynb](https://github.com/KMAnker/IAV_ancestral_reconstruction/blob/main/scripts/h1_xgboost.ipynb) scripts for each protein.
- [figures_tables](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/figures_tables): Figures and tables from the manuscript along with R code used to generate them.
