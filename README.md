# Blueberry rhizosphere fungal communities

Processed data and R scripts associated with the manuscript:

**Host cultivar and geographic region jointly structure blueberry rhizosphere fungal communities with no detectable sequence-variant core**

The dataset comprises 40 ITS1 amplicon samples from eight groups, with five biological replicates per group. Host and site comparisons reuse shared Hefei observations; these observations are not additional replicates and their counts are not summed across contexts.

## Contents

- `inputs/`: processed ASV counts, taxonomy, representative sequences, diversity and distance exports, read-retention records, predicted trophic counts and soil measurements.
- `inputs/original_contexts/`: original host/site comparison inputs and LEfSe exports.
- `tables/`: 36 machine-readable files supporting supplementary Tables S1-S11.
- `code/`: two R scripts and the recorded execution environment.
- `FILE_INDEX.tsv`: supplementary table row counts and fields.
- `README.txt`: detailed statistical conventions, input descriptions and limitations.
- `CHECKSUMS.txt`: SHA-256 checksums; the checksum file itself is excluded.

## Reproduce supplementary statistics

Install R and the packages `vegan`, `DESeq2` (Bioconductor) and `readxl` separately. The recorded environment is in `code/sessionInfo.txt`.

From the repository root:

```sh
Rscript --vanilla code/run_core_diversity.R rerun_core_results
Rscript --vanilla code/run_analysis.R rerun_analysis_results
```

Each script locates the supplied inputs relative to its own path and writes to a new output directory. Supplied manuscript tables are retained as source results. Seeded permutation runs and numerical fitting may produce small differences; consult `README.txt` before comparing rerun values.

## Scope

These scripts reproduce the supplementary statistical analyses specified in `README.txt`. They do not process raw reads or reproduce every manuscript figure. In particular, the combined SIMPER/LEfSe row selection and network role assignment in manuscript Figures 4 and 5 require their original figure-generation sources. Predicted trophic modes are not measurements of biological function. Soil source labels do not establish one-to-one matching to sequencing samples; the supplied soil workbook contains no pH measurements.

## Raw sequence data

BioProject: [PRJNA1551003](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1551003).

The supplied records report SRA submission tracking ID SUB16589303. Formal BioSample/Run accessions and reviewer access details remain pending in the supplied metadata; this repository update does not verify a change to that status. Representative ASV sequences are processed outputs and do not replace raw-read deposition in SRA.

## Citation

Qian H, Feng X, Abdujilil N, Jin J, Cao F, Shi X, Zhou M, Li Y, Zhao T, Zhang Y. 2026. Processed data and R scripts for Host cultivar and geographic region jointly structure blueberry rhizosphere fungal communities with no detectable sequence-variant core. GitHub. https://github.com/Jinxqian-1/blueberry-rhizosphere-fungi.

Use the exact commit or archived version cited by the manuscript when reproducing results. No dataset DOI has been assigned by this update.
