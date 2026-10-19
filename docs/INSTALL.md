# Installation and environment

The workflow combines R packages with independent command-line tools. Install each external tool from its official source and record the exact release or commit in the analysis archive. Do not copy executables or reference panels into Git unless their licenses explicitly permit redistribution.

## Base environment

Create the repository environment:

```bash
conda env create -f environment.yml
conda activate ra-mg-shared-genetics
```

The environment provides R, common R packages, Python, and PLINK. Several research tools use legacy Python 2.7 and should be installed in isolated environments or official containers rather than mixed into the main environment.

## R packages installed outside conda

Use the official installation route for packages that are not reliably available from conda:

```r
install.packages("TwoSampleMR", repos = c(
  "https://mrcieu.r-universe.dev",
  "https://cloud.r-project.org"
))

install.packages("remotes")
remotes::install_github("jean997/cause@v1.2.0")
remotes::install_github("neurogenomics/MAGMA_Celltyping")
```

Install `MRPRESSO` only if MR-PRESSO analysis is required. The MR script saves a clear record when the package is unavailable or the model cannot be fit.

## External tools

Install the tools required by the selected modules and edit `config/config.yml`:

- LDSC and LDSC-SEG
- GNOVA
- ρ-HESS
- GWAS-PW
- MiXeR
- MTAG
- CPASSOC `FunctionSet.R`
- FUSION TWAS
- SMR

GNOVA, ρ-HESS, and MTAG were originally distributed for Python 2.7. Use isolated environments that match the selected release. MiXeR currently recommends an official container for x86 systems; the configuration also supports a source checkout matching the manuscript analysis.

## Reference data

Reference assets are not downloaded automatically. Obtain the population-matched GRCh37 resources specified by each method:

- European HapMap3 SNP list and precomputed LDSC LD scores
- LDSC baseline, regression weights, and tissue `.ldcts` annotations
- European 1000 Genomes PLINK panels for clumping, GNOVA, HESS, MiXeR, FUSION, and SMR
- European LD blocks for HESS and GWAS-PW
- MiXeR LD and random-pruning files
- FUSION GTEx v8 weights and position files
- GTEx v8 and eQTLGen BESD resources for SMR

Confirm that all GWAS coordinates and reference files use GRCh37. Mixing builds can silently remove most variants.

## Local configuration

Create local files that are ignored by Git:

```bash
cp config/config.example.yml config/config.yml
cp config/fusion_weights.example.tsv config/fusion_weights.tsv
cp config/smr_eqtl.example.tsv config/smr_eqtl.tsv
cp config/cpassoc_pruned_snps.example.txt config/cpassoc_pruned_snps.txt
```

Review the cell annotation templates after marker inspection. They are examples from the analysis record, not universal PBMC labels.

## Version capture

Before running the final workflow, save:

```bash
Rscript -e 'writeLines(capture.output(sessionInfo()), "results/R-session-info.txt")'
conda env export --no-builds > results/conda-environment.yml
```

For each external tool, also record `--version` output or the Git commit hash in `results/software-versions.txt`.
