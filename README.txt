SUPPLEMENTAL FILE 2
Host cultivar and geographic region jointly structure blueberry rhizosphere fungal communities with no detectable sequence-variant core

Dataset
Forty ITS1 amplicon samples from eight groups, with five biological replicates per group. The Hefei host comparison contains HFBs, HFNHB, HFRB and HFSHB (20 samples). The site comparison contains HFBs, HFSHB, HNBs, HNSHB, XCBs and XCSHB (30 samples). The shared Hefei samples are reused observations. They are not additional replicates, and their counts are never summed across contexts. The SHB-only comparison contains 15 samples from the same cultivar at the three sites.

Directory structure
inputs: the unified 40-sample ASV counts, taxonomy, representative sequences, alpha-diversity and distance exports, DADA2 read-retention records, aggregated trophic-mode counts, and soil nutrient workbook. original_contexts contains the separate host and site count/taxonomy inputs and their LEfSe marker exports. Hash and ASV-alias feature identifiers are retained because the original ordination and DESeq2 workflows used those respective formats. taxonomic_contrast_plan.csv specifies the fixed rank, taxon and comparison families used for Table S6; it contains no historical test results.
tables: machine-readable values supporting Tables S1-S11 and the supplementary figures. FILE_INDEX.tsv lists file names, row counts and fields. The original input files are not modified by the scripts.
code: two R scripts and the execution environment information. These scripts regenerate the indicated supplementary statistics; they do not redraw or overwrite manuscript figures.

Table conventions
S1: one metadata/QC row per sequencing sample, strict sample-level cores, and pooled-unit ASV intersections. DADA2 read stages and rarefied depth are separate fields. A pooled-unit intersection means presence in at least one sample of every unit; a strict core requires presence in every individual sample.
S2: per-sample diversity indices, group summaries, omnibus Kruskal-Wallis tests, all-pair Dunn tests, and planned Wilcoxon tests for Chao1, Shannon and Simpson diversity (1-D). Dunn and planned tests have separate BH families within each metric and background. Q is the BH-adjusted probability; the redundant unadjusted source P.adj column is omitted.
S3/S4: Hellinger-transformed host/site counts followed by Bray-Curtis distances underpin the main ordinations and matched dispersion tests. The SHB-only comparison uses the supplied unified distance matrix. Original pairwise PERMANOVA exports are retained with SOURCE in their names. Newly seeded permutations are exported separately by the scripts; probabilities from different runs must not be combined in one correction family. The supplied unified-distance host/site diagnostic checks are not part of this supplement.
S5: complete fitted results for the six reported bulk-versus-rhizosphere DESeq2 contrasts, including non-significant rows and missing probabilities. Results_alpha = 0.1 is the independent-filtering optimization parameter. Selected means padj < 0.05 and absolute log2FoldChange > 1. A positive log2FoldChange denotes higher abundance in rhizosphere. The selected counts are 58, 76 and 40 for Hefei NHB, RB and SHB, and 40, 77 and 61 for SHB at Hefei, Huaining and Xuancheng. Prefilter records include all aggregated genera; complete fitted results include only those passing count-sum > number of samples. The additional non-bulk host contrasts and alpha = 0.05 sensitivity analyses are excluded.
S6: genus/phylum per-sample proportions are fractions; test-table means and SD are percentages. Full ITS libraries, including unclassified reads, form the denominators. BH adjustment is separate by rank, background and question. Family_n counts valid tests; an all-zero untestable contrast retains a missing P/Q value. FoldChange = mean(group 2)/mean(group 1), and is missing if the first mean is zero.
S7: predicted trophic-mode percentages, means and planned tests. Mixed annotations are assigned in the order pathotroph, symbiotroph, saprotroph. Unassigned reads remain in library denominators. Planned tests use nine hypotheses per background. S7b uses three omnibus mode tests per compartment; S7c uses nine pairwise site-by-mode tests per compartment. These are different BH families. Predictions are not measured functions, and annotation-confidence categories are unavailable in the supplied aggregated output.
S8: site descriptions and coordinates.
S9/S10: the 40 measured soil C/N/P/K source records, summaries and coefficients of variation. Units are g kg-1; SD uses n-1 and CV = 100 x SD/mean. Original workbook labels, sheets and row numbers are retained. These labels do not establish one-to-one matching to sequencing samples. The workbook contains no pH values. Hefei rhizosphere summaries pool 15 records across three hosts; the host-level table retains five records per group.
S11: full host/site LEfSe marker exports and the LDA >= 4 display subset. Identical P.unadj and P.adj source columns are reported as nominal probabilities, not established multiple-testing-adjusted results. Marker-only exports do not represent all tested taxa.
Missing numeric cells are missing values, not zero or evidence of no effect. Source identifiers and taxonomic labels are retained.

Running the scripts
Dependencies: R, vegan, DESeq2 and readxl. Tested with R 4.6.0, vegan 2.7-5, DESeq2 1.52.0 and readxl 1.5.0; full environment information is in code/sessionInfo.txt. Install dependencies separately; the scripts do not install packages or download data.
From the extracted Supplemental_data directory:
  Rscript --vanilla code/run_core_diversity.R /path/to/new_core_results
  Rscript --vanilla code/run_analysis.R /path/to/new_analysis_results
On Windows, use an available UTF-8 locale; the scripts set LC_CTYPE when possible. Both scripts locate inputs relative to their own files. Output defaults to rerun_results; supply separate output directories to keep executions distinct. The first script exports alpha-diversity, ASV-core/intersection and SHB-only community checks. The second exports taxonomic tests, Hellinger community/dispersion results, predicted trophic-mode tests, the six original-parameter DESeq2 fits, soil summaries and LEfSe source subsets. The supplied tables retain the manuscript source results. New permutation runs can produce small Monte Carlo differences in P values; fixed effect statistics and interpretation should be checked separately. Minor numerical optimizer differences in DESeq2 do not change the selected-genus sets in the verified execution.

Scope
The two scripts reproduce supplementary statistical tables from the supplied processed inputs. They are not raw-read processing pipelines and do not recreate the combined SIMPER/LEfSe row selection or network role assignment in manuscript Fig. 4 and 5. Those figure workflows require their corresponding original figure-generation sources. This data supplement does not add experiments, independent samples, measured functional data or pH records.

Sequence archive
BioProject PRJNA1551003; SRA submission tracking SUB16589303. Formal BioSample/Run accessions and reviewer access details remain pending, as stated in the manuscript. Representative ASV sequences in this supplement do not replace deposition of raw reads in SRA.

Repository
https://github.com/Jinxqian-1/blueberry-rhizosphere-fungi
The repository preserves the inputs/, tables/ and code/ structure used by both scripts.
