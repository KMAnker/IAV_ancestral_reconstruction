
# Exploring genetic signatures of zoonotic influenza A virus at the swine-human interface with phylogenetic and ancestral sequence reconstruction

Pipeline, scripts, and data for the paper:
Anker, K.M., Ciucani, M.M., Nissen, J.N., Anderson, T.K., Pedersen, A.G., and Trebbien, R. (in preparation). Exploring genetic signatures of zoonotic influenza A virus at the swine-human interface with phylogenetic and ancestral sequence reconstruction.


## Analysis pipeline and scripts:
The steps of the pipeline for data preparation and analysis is described in the [ancestor_project_notes](https://github.com/KMAnker/IAV_ancestral_reconstruction/blob/main/ancestor_project_notes.sh) file.
Here, it can also be seen where each of the scripts in the [scripts](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/scripts) folder is used


## Data:
Sequence data for this analysis was obtained from GISAID and NCBI Influenza Virus Resource (Genbank). 
We gratefully acknowledge all data contributors, i.e., the authors and their originating laboratories responsible for obtaining the specimens, and their submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID initiative as well as the NCBI Influenza Virus Resource database, on which this research is based.

We cannot publicly upload the genetic sequences from GISAID, and therefore no sequence (fasta) files are uploaded to this repository. However, the [data](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/data) folder contains lists of accession numbers for all the sequences originally downloaded for analysis.


## Analysis and output folders:
Folders containing input and output files for the different steps of the analysis:
- [treefiles](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/treefiles): IQtree output files and both input and output files from Treetime analysis (except fasta sequence files)
- [anclib_files](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/anclib_files): Output files from the anclib analysis
- [selection_analysis](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/selection_analysis): Output json files from the HyPhy aBSREL and MEME analyses as well as text and excel files containing concatenated output from all segments
- [xgboost](https://github.com/KMAnker/IAV_ancestral_reconstruction/tree/main/xgboost): one-hot encoded input files and output files from xgboost analyses (run with [xgboost.ipynb](https://github.com/KMAnker/IAV_ancestral_reconstruction/blob/main/scripts/h1_xgboost.ipynb) scripts.
 

## Small benchmark 


<sup>1</sup> Criteria for convergence for simulated Sq data is defined as the number of steps until Rwp < 0.04.<br>
<sup>2</sup> Criteria for convergence for simulated Gr data is defined as the number of steps until Rwp < 0.04.<br>
<sup>3</sup> Criteria for convergence for experimental Sq data is defined as the number of steps until Rwp < 0.84.<br>
<sup>4</sup> Criteria for convergence for experimental Gr data is defined as the number of steps until Rwp < 0.80.<br>
<sup>5</sup> Criteria for convergence for multi-objective optimisation is defined for both above criteria.

# Usage
See (https://github.com/AndySAnker/ScattBO/tree/main/tools) for examples of single-objective optimisation with [Dragonfly](https://github.com/dragonfly/dragonfly/tree/master) or [skopt](https://scikit-optimize.github.io/stable/auto_examples/bayesian-optimization.html).

## Example usage with [Dragonfly](https://github.com/dragonfly/dragonfly/tree/master)


## Visualise results
```python

# Use the parameters that minimize the function as input to the ScatterBO_small_benchmark function
ScatterBO_small_benchmark((pH, pressure, solvent), plot=True, simulated_or_experimental='simulated', scatteringfunction='Gr')

```


## Reporting issues

If you encounter any issues or problems with our software, please report them by opening an issue on our GitHub repository. Please include as much detail as possible, including steps to reproduce the issue and any error messages you received.

## Seeking support

If you need help using our software, please reach out to us on our GitHub repository. We'll do our best to assist you and answer any questions you have.

# References

See also the [Twitter submission post](https://twitter.com/SodeAndy/status/1773474538631651769)
