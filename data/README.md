# Data retrieval and local filenames

Large public matrices are intentionally excluded from version control. Download them from the cited repositories and store them outside Git, or under a local ignored `data/raw/` directory.

The legacy scripts expect the following working files:

| Expected filename/path | Description |
|---|---|
| `/tmp/gtex_median.gct.gz` | GTEx v8 tissue-median gene TPM matrix |
| `/tmp/gtex_gene_tpm.gct.gz` | GTEx v8 specimen-level gene TPM matrix |
| `/tmp/human_meth.txt.gz` | GSE213478 processed methylation matrix |
| `/tmp/human_series.txt.gz` | GSE213478 GEO series metadata |
| `/tmp/epic_manifest.csv.gz` | EPIC hg38 probe manifest |
| `/tmp/mouse_bodymap.xls` | Li et al. mouse BodyMap supplementary expression workbook |
| `/tmp/mouse_meth.csv.gz` | GSE290585 processed mouse methylation matrix |
| `/tmp/mouse_series.txt.gz` | GSE290585 GEO series metadata |
| `/tmp/mm10.gtf.gz` | Mouse gene annotation used for region mapping |
| `/tmp/ensembl_human_mouse.tsv` | Ensembl one-to-one orthologue export |
| GSE92486 RNA-seq/WGBS files | Independent paired mouse-liver replication inputs |

Do not commit controlled-access donor covariates, credentials, or raw GTEx individual-level files. Verify checksums and source versions after download, and record them in a release manifest.

