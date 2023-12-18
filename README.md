# Scripts for computational analyses of scRNA-seq and single-cell multiome data presented in "Partitioning heritability using single-cell multi-omics identifies a novel macrophage subpopulation conveying genetic risks of coronary artery disease"

## Description

This repository contains the scripts to reproduce results and figures for the multiome paper:

Jiang J, Hiron T.K., et al. "Partitioning heritability using single-cell multi-omics identifies a novel macrophage subpopulation conveying genetic risks of coronary artery disease." bioRxiv (2023): 2023-09.
https://www.biorxiv.org/content/10.1101/2023.09.14.557845v1

## Table of Contents

- [Dependencies](#dependencies)
- [Content](#content)
- [License](#license)

## Dependencies
To run the scripts in this repository, first make sure that all packages listed in the [DEPENDS_R](DEPENDS_R.txt) file have been installed. 
Additional utility functions can be found here(https://github.com/jhjiang2020/utils_R).
To replicate the deep-learning based macrophage classification results, you may need to install the Python dependencies as listed in [scvi-tools](https://github.com/scverse/scvi-tools).

## Content
Here is a breakdown of all scripts in this repository by dataset and analysis:
1. The `scRNA-seq` folder contains scripts for sample integration and clustering for *ex vivo* macrophages(scRNA-seq/001_Seurat_CCA_Integration.Rmd) and for pseudo-bulk DE analysis(scRNA-seq/002_Pseudobulk_DE.Rmd).
2. The `multiome` folder contains all scripts from integration(multiome/001_WNN_integration.Rmd), DE(multiome/002_multiome_DE.Rmd), CRE mapping(multiome/003_CRE_mapping.Rmd), GWAS variants prioritization(multiome/004_GWAS_SNPs.Rmd), and s-LDSC analysis(multiome/005_s-LDSC_analysis.Rmd).
3. The `meta-analysis` folder contains scripts for sample integration and clusterings at three different levels (all cells(meta-analysis/001_meta_human_plaques.Rmd), macrophages(meta-analysis/002_meta_macrophages.Rmd) and LAMs(meta-analysis/004_cluster_stability_LAM.Rmd)). We also include Python scripts for the deep-learning-based macrophage classification(meta-analysis/003_scVI_labeltransfer.ipynb).

**A special note for running s-LDSC analysis:**
We performed our s-LDSC analysis as instructed in the [ldsc wiki](https://github.com/bulik/ldsc/wiki), using the reference 1000G phase 3 genome (EUR) in hg38 obtained from [PLINK resource page](https://www.cog-genomics.org/plink/2.0/resources). The script we used for running s-LDSC was not included in this repository as it was specifically designed for our high-performance biomedical computing cluster. The LDSC code we executed is outlined below:
```bash
conda activate ldsc-env

python ldsc/make_annot.py \
--bed-file $BEDFile \
--bimfile $Reference \
--annot-file $OutDir

python ldsc/ldsc.py \
--l2 \
--bfile $Reference \
--ld-wind-kb 1000 \
--annot $AnnotFile \
--thin-annot \
--out ${OutDir}

python ldsc/ldsc.py \
--h2 $sumstats \
--overlap-annot  \
--frqfile $Reference \
--w-ld $RefWeight \
--ref-ld $Annot1,$Annot2 \
--out $OutName
```


## License

This repository is licensed under the [MIT License](LICENSE).