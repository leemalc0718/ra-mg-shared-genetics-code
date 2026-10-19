# Shared genetic architecture of rheumatoid arthritis and myasthenia gravis

This repository contains analysis code accompanying a study of the shared genetic architecture of rheumatoid arthritis (RA) and myasthenia gravis (MG). It organizes the manuscript workflow into auditable modules for global and local genetic correlation, polygenic overlap, causal inference, tissue and cell-type enrichment, cross-trait association, colocalization, gene prioritization, and single-cell expression mapping.

No individual-level data, licensed reference panels, GWAS summary statistics, credentials, analysis results, or institution-specific filesystem paths are included. All paths are supplied through `config/config.yml`, created locally from `config/config.example.yml`.

## Repository structure

```text
config/       Example configuration, loci, and annotation tables
data/         Local input area; ignored by Git
docs/         Installation, input schema, method mapping, and official resources
R/            Shared R utilities
scripts/      Numbered analysis modules
tests/        Static release checks
results/      Generated outputs; ignored by Git
```

## Quick start

1. Clone the repository and copy the example configuration.

   ```bash
   cp config/config.example.yml config/config.yml
   ```

2. Install the R environment and the external tools listed in [docs/INSTALL.md](docs/INSTALL.md). Record exact releases or commit hashes before analysis.

3. Convert each GWAS to the standardized schema.

   ```bash
   Rscript scripts/00_standardize_gwas.R \
     --input data/raw/ra.tsv.gz \
     --output data/processed/ra.tsv.gz \
     --snp SNP --chr CHR --bp BP --a1 A1 --a2 A2 \
     --beta beta --se se --p p --eaf freq --n n --info INFO
   ```

4. Update `config/config.yml`, then validate the inputs.

   ```bash
   Rscript scripts/00_validate_inputs.R config/config.yml
   ```

5. Run only the modules needed for the analysis. For example:

   ```bash
   bash scripts/01_ldsc.sh config/config.yml
   Rscript scripts/06_two_sample_mr.R config/config.yml ra mg
   Rscript scripts/12_coloc.R config/config.yml config/loci.tsv
   ```

The scripts do not download controlled or large reference datasets automatically. This prevents accidental redistribution and makes provenance explicit.

## Manuscript-specific analysis map

| Analysis | Module |
|---|---|
| Global genetic correlation | `01_ldsc.sh`, `02_gnova.sh` |
| Local genetic correlation | `03_hess.sh`, `04_gwas_pw.sh` |
| Polygenic overlap | `05_mixer.sh` |
| Bidirectional MR and CAUSE | `06_two_sample_mr.R`, `07_cause.R` |
| Tissue and cell-type enrichment | `08_ldsc_seg.sh`, `09_magma_celltyping.R` |
| Cross-trait association | `10_mtag.sh`, `11_cpassoc.R` |
| Bayesian colocalization | `12_coloc.R` |
| TWAS and SMR | `13_fusion_twas.sh`, `14_smr.sh` |
| Single-cell expression mapping | `15_scrna_ra.R`, `16_scrna_mg.R` |

Detailed thresholds and method-to-code correspondence are in [docs/METHODS_TO_CODE.md](docs/METHODS_TO_CODE.md).

## Reproducibility notes

- The standardized GWAS schema uses GRCh37 coordinates and effect estimates aligned to `A1`.
- The original analysis excluded the MHC region for MiXeR and used European ancestry LD references.
- `coloc.abf` receives all harmonized variants in each region; regional variants are not clumped.
- Cluster-to-cell-type mappings are stored in reviewable TSV files. Scripts stop when observed cluster IDs do not match the supplied mapping.
- Manual deletion of SNP rows is not reproducible and is not used. MR-PRESSO outliers, when available, are identified from model output and recorded.

## Validation

Run the release checks before publishing:

```bash
bash tests/static_checks.sh
```

These checks parse all R scripts, validate shell syntax and YAML, and reject personal paths, temporary R paths, and credential-like strings.

## Data availability

GWAS, eQTL, single-cell, and reference-panel data must be obtained from the original studies or official repositories under their applicable terms. See [docs/DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) for the required local formats.

## Citation

When using this code, cite the associated article and each underlying method or dataset. Official software resources are listed in [docs/OFFICIAL_RESOURCES.md](docs/OFFICIAL_RESOURCES.md).
