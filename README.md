# ARPHI-Laguna-Lake-Stat-Analysis

This repository contains the R scripts and analytical workflows used to process and statistically analyze antibiotic resistance gene (ARG) occurrence in environmental *Escherichia coli* isolates and physicochemical water-quality data from Laguna Lake and its tributary rivers, Philippines.

## Data Standardization and Physicochemical Data Processing

The `raw_data/` folder contains the raw encoded data for bacterial isolates and physicochemical records obtained from the Laguna Lake Development Authority (LLDA).

Data standardization was performed using [`standardize_data.R`](codes/standardize_data.R). This script generates the standardized datasets used in the subsequent analyses:

- **[`ARPHI-isolate-database.tsv`](ARPHI-isolate-database.tsv)** — Contains all presumptive *E. coli* isolates collected throughout the study, including uidA confirmation results and ARG PCR results.
- **[`sample_metadata.tsv`](sample_metadata.tsv)** — Contains metadata for each independent water sample, including collection date, sampling site, total number of isolates, number of isolates retained for analysis, and physicochemical parameters.

## ARG Occurrence, Co-occurrence, and Pairwise Associations

[`ARG_association_analysis.R`](codes/ARG_association_analysis.R) generates the ARG occurrence and co-occurrence visualizations, including the UpSet plot (Figure 2). It also performs pairwise ARG association analyses and generates the corresponding odds ratio forest plot (Figure 3).

## Spatial and Seasonal Patterns of ARG Prevalence

[`ARG_pattern.R`](codes/ARG_pattern.R) evaluates the associations of geographic region and season with ARG prevalence using generalized linear mixed-effects models (GLMMs) with a binomial error distribution. The script also generates the estimated marginal means (EMMs) plot (Figure 4).

## Physicochemical Data Analysis

[`pc_analysis.R`](codes/pc_analysis.R) processes the physicochemical dataset by substituting below-reporting-limit values with the reporting limit divided by √2 and imputing missing measurements using multiple imputation by chained equations (MICE). It subsequently performs principal component analysis (PCA) and multivariate analyses, including linear mixed-effects modeling (LMM) of PC1 scores and PERMANOVA. This script also generates the PCA plot (Figure 5).

## Associations Between ARG Prevalence and Physicochemical Parameters

[`ARG_PC_correlation.R`](codes/ARG_PC_correlation.R) evaluates associations between ARG prevalence and physicochemical parameters using beta-binomial generalized linear mixed-effects models. The script also generates the odds ratio forest plot (Figure 6).
