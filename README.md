# Tissue specificity and methylation–expression coupling

Reproducibility package for the comparative human–mouse analysis of gene-expression tissue specificity and promoter/gene-body DNA methylation–expression coupling.

## Study design

The analysis has two complementary levels:

1. **Tissue-average (ecological) analysis:** gene-level methylation and expression are correlated across harmonized tissue means.
2. **Same-biospecimen validation:** methylation and RNA measurements from the same specimens are analyzed within human colon, lung and skeletal muscle and independently in mouse liver.

Tissue specificity is quantified with the tau index. The principal fixed groups are `tau < 0.4` (broad expression) and `tau > 0.8` (tissue-restricted expression). Human–mouse comparisons use Ensembl one-to-one orthologues and a matched 15-organ expression panel.

## Repository layout

```text
config/          Dataset accession and analysis settings
data/            Retrieval and local-file instructions (no large raw data)
docs/            Analysis notes and manuscript-facing documentation
manuscript/      Journal manuscript and editable figure deck
results/         Derived workbooks suitable for audit and reuse
scripts/         Python analysis and figure-generation scripts
submission/      Elsevier highlights and pre-submission checklist
```

## Public datasets

| Component | Species | Resource |
|---|---|---|
| Expression atlas | Human | GTEx v8 median gene TPM |
| DNA methylation atlas and paired validation | Human | GEO `GSE213478` |
| Expression atlas | Mouse | Li et al. 17-tissue BodyMap |
| DNA methylation atlas | Mouse | GEO `GSE290585` |
| Independent paired RNA-seq/WGBS replication | Mouse | GEO `GSE92486` |
| One-to-one orthology | Human/mouse | Ensembl Compara |

See `data/README.md` for the expected local filenames. Source data are not redistributed.

## Environment

Create the analysis environment with either:

```bash
conda env create -f environment.yml
conda activate tsi-methylation
```

or:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Reproduction notes

The archived scripts reproduce the exact analysis history and retain some legacy `/tmp` paths used during the original run. Before a public release, set up the filenames listed in `data/README.md` or replace those path constants with a project-specific data directory. Fixed random seeds are documented in the manuscript and scripts. Derived workbooks are included so every reported summary can be audited without downloading the large source matrices.

## Principal outputs

- `results/Human_Mouse_TSI_Methylation_Correlation.xlsx`
- `results/Human_Mouse_TSI_Tissue_Overlap.xlsx`
- `results/Human_Mouse_TSI_Validation_Analyses.xlsx`
- `manuscript/Tissue_Specificity_Methylation_Expression_CBP_D_Submission.docx`
- `manuscript/TSI_Manuscript_Figures_Editable.pptx`

## Citation

Use `CITATION.cff` after completing the author list, repository URL, release date and DOI. Cite the original public datasets separately as listed in the manuscript.

## License

No reuse license has been assigned in this pre-submission package. The authors should select an explicit code and data license before making the repository public.

