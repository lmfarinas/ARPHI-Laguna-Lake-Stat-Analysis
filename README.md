# ARPHI-Laguna-Lake-Stat-Analysis

## Data Standardization and Physicochemical Data Processing
The raw_data/ folder contains the raw encoded data for bacterial isolates, and encoded LLDA physicochemical records.

Data standardization was performed using [`standardize_data.R`](codes/standardize_data.R). This produces the data used for the analyses:
- **[`ARPHI-isolate-database.tsv`](ARPHI-isolate-database.tsv)** — Contains all presumptive *E. coli* isolates collected throughout the study, including uidA confirmation and antibiotic resistance gene (ARG) PCR results.
- **[`sample_metadata.tsv`](sample_metadata.tsv)** — Contains metadata for each independent water sample, including collection date, sampling site, total number of isolates, number of isolates retained for analysis, and physicochemical parameters.
