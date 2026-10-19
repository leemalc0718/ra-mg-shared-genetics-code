# Data dictionary

## Standardized GWAS summary statistics

All downstream R modules expect a tab-separated plain-text or gzip-compressed file with one row per unique variant.

| Column | Type | Definition |
|---|---|---|
| `SNP` | character | dbSNP rs identifier |
| `CHR` | integer | Autosomal chromosome, 1 through 22 |
| `BP` | integer | One-based GRCh37 base-pair position |
| `A1` | character | Effect allele for `beta` and `Z` |
| `A2` | character | Other allele |
| `beta` | numeric | Log-odds or linear-regression effect estimate |
| `se` | numeric | Standard error of `beta` |
| `p` | numeric | Two-sided association P value |
| `eaf` | numeric | Effect-allele frequency in `(0,1)` |
| `n` | numeric | Per-variant or constant total sample size |
| `info` | numeric | Imputation INFO score; optional at import, required for the documented INFO filter |
| `Z` | numeric | Signed Z score computed as `beta / se` |

`00_standardize_gwas.R` maps source column names to this schema, uppercases alleles, computes `Z = beta / se`, applies configured quality filters, and rejects duplicate SNP identifiers.

## Case-control metadata

Case and control counts are stored under `traits` in `config/config.yml`. The observed case fraction passed to `coloc` is calculated as:

```text
s = n_case / (n_case + n_control)
```

Population prevalence is not inferred. Add an external prevalence estimate only when a selected method explicitly requires liability-scale conversion.

## Method-specific derived inputs

- **LDSC:** generated with the official `munge_sumstats.py` and restricted to HapMap3 alleles.
- **GNOVA:** use the format required by the installed GNOVA release.
- **ρ-HESS:** `SNP CHR BP A1 A2 Z N`.
- **GWAS-PW:** sorted `SNPID CHR POS Z_RA V_RA Z_MG V_MG`, with alleles harmonized before construction. Here `V_*` is the variance of the effect estimate.
- **MiXeR:** LDSC-compatible `.sumstats.gz` input.
- **FUSION:** whitespace-delimited summary statistics containing at least `SNP`, `A1`, `A2`, and `Z`.
- **SMR:** GWAS input formatted according to the SMR documentation and eQTL data supplied as a matching BESD prefix.

## Single-cell inputs

### RA

Each sample directory must contain a standard 10X matrix triplet:

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

The first feature column is treated as the gene identifier. If gene symbols are in the second column, change the `gene.column` setting in the script before analysis and record the decision.

### MG

Each input is a tab-separated gene-by-cell count matrix. The first column contains gene names and all remaining columns contain non-negative counts. File names determine the sample and group using the configured disease and control patterns.

## Data that must not be committed

The `.gitignore` excludes raw and processed GWAS files, reference panels, single-cell matrices, generated results, logs, R objects, and local configuration. Before publishing, run `tests/static_checks.sh` and inspect `git status`.
