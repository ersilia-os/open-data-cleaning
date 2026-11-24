# Open data cleaning

[![status](https://img.shields.io/badge/status-ready-green)](https://github.com/)

Notebooks that clean and curate open datasets downloaded from public sources.

## Background

Ersilia collects public bioactivity and compound-library data to train and benchmark models. This
repository keeps the original downloads, the cleaning notebooks and the cleaned datasets, so each
dataset can be traced back to its source. All data here is public. Data shared by collaborators
does not belong in this repository.

## Datasets

Original downloads are in `data/raw/` and cleaned files in `data/processed/`, one folder per dataset.

| Dataset | Content | Source | Notebook |
|---|---|---|---|
| `cyp450` | CYP2C9, CYP2C19, CYP2D6 and CYP3A4 inhibition, merged into `processed/cyp450/cyp_all_<cyp>.csv` | PubChem assays 1851, 410, 883, 884, 891 and 899; ChEMBL IC50 and Ki | `cyp450.ipynb` |
| `esaki2019` | Metabolic stability (stable or unstable) of ChEMBL compounds | Esaki *et al.*, *Molecular Informatics* 2019 | `esaki2019.ipynb` |
| `mmv_boxes` | Compound lists of the Pathogen Box (with *M. tuberculosis*, *P. falciparum* and HepG2 data), Malaria Box, Pandemic Response Box and Global Health Priority Box | Medicines for Malaria Venture (MMV) | `mmv_boxes.ipynb` |
| `kinase_prom.csv` | SMILES with a binary activity label | Unknown: no notebook or record of origin | none |

The final CYP files (`cyp_all_<cyp>.csv`) were used to train the CYP models of the
[H3D screening cascade](https://github.com/ersilia-os/h3d-screening-cascade-models) (eos7kpb).

**Known gaps:** some notebooks were edited after the stored files were made, so a few file names
do not match. `cyp450.ipynb` writes `pr_aid1851_*`, `pr_pubchem_*` and `pr_cyp_all_*`, while the
stored files are `pr_AID_*`, `pr_PUBCHEM*` and `cyp_all_*`. `esaki2019.ipynb` writes
`pr_esaki2021.csv`, while the stored file is `pr_esaki2019.csv`. The stored files are the reference.

## History

Until November 2026 this repository stored its data with Git LFS. It was recreated from the
analysis template without Git LFS. Every file is well under GitHub's 100 MB limit, so the data is
kept in git. The CO-ADD and quorum-sensing (QSI) biofilm datasets
were dropped at that point.

## Repository structure

The folder layout is described in [`CLAUDE.md`](CLAUDE.md). Git tracks the notebooks and the data, as no file is too large for GitHub. Git LFS is not used.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
