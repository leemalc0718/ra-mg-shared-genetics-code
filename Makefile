.PHONY: check validate ldsc gnova hess gwas-pw mixer mr cause ldsc-seg magma mtag cpassoc coloc twas smr scrna-ra scrna-mg genes

CONFIG ?= config/config.yml

check:
	bash tests/static_checks.sh

validate:
	Rscript scripts/00_validate_inputs.R $(CONFIG)

ldsc:
	bash scripts/01_ldsc.sh $(CONFIG)

gnova:
	bash scripts/02_gnova.sh $(CONFIG)

hess:
	bash scripts/03_hess.sh $(CONFIG)

gwas-pw:
	bash scripts/04_gwas_pw.sh $(CONFIG)

mixer:
	bash scripts/05_mixer.sh $(CONFIG)

mr:
	Rscript scripts/06_two_sample_mr.R $(CONFIG)

cause:
	Rscript scripts/07_cause.R $(CONFIG)

ldsc-seg:
	bash scripts/08_ldsc_seg.sh $(CONFIG)

magma:
	Rscript scripts/09_magma_celltyping.R $(CONFIG)

mtag:
	bash scripts/10_mtag.sh $(CONFIG)

cpassoc:
	Rscript scripts/11_cpassoc.R $(CONFIG)

coloc:
	Rscript scripts/12_coloc.R $(CONFIG) config/loci.tsv

twas:
	bash scripts/13_fusion_twas.sh $(CONFIG)

smr:
	bash scripts/14_smr.sh $(CONFIG)

scrna-ra:
	Rscript scripts/15_scrna_ra.R $(CONFIG)

scrna-mg:
	Rscript scripts/16_scrna_mg.R $(CONFIG)

genes:
	Rscript scripts/17_prioritize_genes.R $(CONFIG)
