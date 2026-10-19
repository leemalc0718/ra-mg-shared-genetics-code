# Methods to code mapping

This document maps the manuscript methods to repository modules and separates manuscript-defined thresholds from software defaults.

| Manuscript component | Code | Main settings |
|---|---|---|
| Summary-statistic QC | `00_standardize_gwas.R` | GRCh37, autosomes, MAF at least 0.01, INFO at least 0.90 when available |
| Global genetic correlation | `01_ldsc.sh`, `02_gnova.sh` | European LD reference; LDSC HapMap3 munging |
| Local genetic correlation | `03_hess.sh` | 22 autosomes; European LD blocks; sample overlap and phenotypic correlation are explicit configuration values |
| Shared local regions | `04_gwas_pw.sh` | Harmonized signed effects; recommended LD blocks; model 3 posterior probability above 0.80 used downstream |
| Polygenic overlap | `05_mixer.sh` | MHC exclusion at chr6:26-34 Mb; 20 random SNP subsets; univariate and bivariate fit and test stages |
| Bidirectional MR | `06_two_sample_mr.R` | P below 5e-8, r2 below 0.001 within 10,000 kb, F at least 10, harmonization action 2 |
| Correlated pleiotropy | `07_cause.R` | Nuisance estimation from at most one million variants; exposure clumping at P below 1e-3 and r2 below 0.01 |
| Tissue enrichment | `08_ldsc_seg.sh` | LDSC-SEG with GTEx multi-tissue annotations; FDR is calculated from output during reporting |
| Cell-type enrichment | `09_magma_celltyping.R` | 35 kb upstream and 10 kb downstream; linear and top-10-percent tests; persistent cache |
| Cross-trait MTAG | `10_mtag.sh`, `10b_mtag_postprocess.R` | MAF at least 0.01; N at least 75 percent of the 90th percentile; genome-wide threshold 5e-8; 500 kb, r2 0.2 clumping; nominal original-trait evidence and concordant directions |
| Cross-trait CPASSOC | `11_cpassoc.R` | SHet statistic; 20 folds; threshold 5e-8; both original traits below 1e-3; 500 kb, r2 0.2 clumping |
| Colocalization | `12_coloc.R` | Dense harmonized variants within 250 kb of nine lead variants; case-control proportions; PP.H4 above 0.80 |
| TWAS | `13_fusion_twas.sh`, `17_prioritize_genes.R` | GTEx v8 weights for whole blood, lung, spleen, and EBV-transformed lymphocytes; within-analysis BH FDR below 0.05 |
| SMR and HEIDI | `14_smr.sh`, `17_prioritize_genes.R` | GTEx v8 and eQTLGen; BH FDR below 0.05; HEIDI P above 0.05; more than 10 HEIDI SNPs when reported |
| RA and MG single-cell expression | `15_scrna_ra.R`, `16_scrna_mg.R` | LogNormalize, 2,000 variable genes, PCA, Louvain resolution 0.5, UMAP; 10 PCs for RA and 15 for MG |

## Deliberate corrections to the analysis notes

- A credential embedded in the original MAGMA notes is omitted. Authentication, if needed for a private resource, must use a local credential manager or environment variable outside Git.
- Temporary R session directories are not reused as stable MAGMA inputs.
- MR-PRESSO outliers are obtained from model output. No SNPs are removed by hard-coded row number or visual inspection.
- CAUSE samples `min(n_variants, 1,000,000)` variants for nuisance estimation.
- MiXeR includes both fit and test stages for univariate and bivariate models.
- `coloc.abf` receives a dense regional set rather than LD-pruned lead variants.
- Single-cell preprocessing is performed once after combining QC-filtered count matrices. Cluster labels are loaded from a versioned table and must match observed cluster IDs.

## Analyses requiring manual scientific review

The repository automates computation but does not replace the following decisions:

- evaluation of GWAS ancestry, build, sample overlap, and phenotype definitions;
- selection of HESS phenotypic-correlation correction when samples overlap;
- inspection of MiXeR model-fit diagnostics;
- review of MR instrument validity and residual pleiotropy;
- review of MAGMA multiple-testing denominators for each cell-type dataset;
- confirmation of MTAG maxFDR assumptions;
- assessment of the single-causal-variant assumption for `coloc.abf`;
- marker-based validation of all single-cell annotations.
