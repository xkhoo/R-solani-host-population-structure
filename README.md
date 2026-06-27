# _Rhizoctonia solani_ host-associated population structure

This repository contains analysis scripts, metadata extraction scripts, output files (mainly from metadata extraction), and representative diagnostic plots associated with a study of host-associated population structure and host-transition patterns in the *Rhizoctonia solani* species complex.

The analyses focus on public *R. solani* sequence records and associated metadata, with particular emphasis on host-associated genetic structure, population differentiation, stratified host-group subsampling, and phylogenetic host-transition analyses.

## Data and metadata

The `data/` directory contains accession information used to retrieve sequence metadata from GenBank.

The `output/` directory contains metadata tables extracted from GenBank records, including host, location, collection date, and strain information. These files are generated using the Python scripts in `scripts/metadata_extraction/`.

The original sequence alignments and large intermediate phylogenetic files are not included in this repository because of file size. Sequence accessions are provided so that raw records can be retrieved from NCBI GenBank, and their associated metadata are included in the Supplementary Table S1 of the manuscript. Multiple sequence alignments and large intermediate files are available from the author upon reasonable request.

## Scripts

### Metadata extraction

The scripts in `scripts/metadata_extraction/` use Biopython Entrez utilities to retrieve metadata from GenBank records based on accession numbers.

These scripts extract:

* host metadata
* geographic location metadata
* collection date metadata
* strain/isolate metadata

Example:

```bash
python scripts/metadata_extraction/extract_host_data_from_gb.py
```

The scripts assume that the accession list is stored in:

```text
data/accessions.txt
```

and write extracted metadata files to:

```text
output/
```

### Population structure analyses

The repository includes R scripts for population-structure and differentiation analyses.

`DAPC_population_structure_Rsolani.R` performs discriminant analysis of principal components using sequence-derived genetic data and associated metadata.

`AMOVA.R` contains distance-based population differentiation analyses, including pairwise comparisons among host or geographic groups.

### Stratified host-group subsampling

`Stratified_sampling_hostgroup.R` generates host-stratified subsamples for downstream Markov-jumps analyses. This was used to reduce host-group sampling imbalance by drawing comparable numbers of sequences from major host families. **10**

The major host families considered in the stratified subsampling workflow include:

```text
Solanaceae
Amaranthaceae
Fabaceae
Brassicaceae
Poaceae
```

### Host-transition network analysis

`Host_transition_network_analysis.R` summarizes host-transition counts from BEAST infered host-transition output files and generates host-transition network visualizations.

Large BEAST input and output files are not stored in this GitHub repository. To reproduce this analysis, users need to generate or obtain the corresponding BEAST annotated-history files before running the script.

## Supplementary diagnostics

The `supplementary diagnostics/` directory contains representative BEAST diagnostic outputs. These include selected Tracer plots and LogAnalyser summaries used to evaluate convergence, effective sample sizes, and posterior similarity between independent chains for representative empirical subsamples.

These files are included as examples of MCMC diagnostic checks for representative empirical subsamples. They are not intended to represent the full set of diagnostic plots generated for every subsample and chain.

## Notes on reproducibility

Some scripts contain project-specific file paths from the original working environment. Users may need to update file paths before running the scripts on a different machine.

This repository is intended to document the analysis workflow and provide reusable scripts and metadata outputs associated with the study. Large sequence alignments, BEAST XML files, tree files, and annotated-history files are excluded from the GitHub repository but can be made available upon reasonable request.

## Citation

If using scripts or metadata from this repository, please cite the associated manuscript once available. For questions about the repository, sequence alignments, or intermediate analysis files, please contact the repository owner.
