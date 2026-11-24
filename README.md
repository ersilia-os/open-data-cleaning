# Ersilia's analysis template

[![status](https://img.shields.io/badge/status-pending-red)](https://github.com/)

This repository provides a structured template for setting up new research analysis in Ersilia.

## Background

Replace this paragraph with a short description of the project. This description should explain the background or context of the project, specifying collaborators.

## Project tracking

* [Meeting #1](https://example.com) - YYYY-MM-DD: Short description
* ...

## Progress

### Step 1

Describe what was done in step 1.

### Step 2

Describe what was done in step 2.

## Repository structure

This repository is organized as follows:

- **data/**
  - **raw/** → Original, untouched datasets  
  - **processed/** → Cleaned and transformed datasets  

- **scripts/** → Standalone scripts for preprocessing or automation. Numbered in sequential order for running.

- **notebooks/** → Jupyter notebooks for exploration and prototyping.

- **assets/** → Images, figures, and other static resources  

- **output/** → Results of the scripts, numbered by file or folder according to the scripts numbering  

- **src/** → Core source code and reusable modules  

- **tools/** → Helper utilities and development tools  

- **docs/** → Project documentation and reports.

- **tmp/** → Temporary files or intermediate outputs  

The project is tracked in [GitHub](https://github.com/ersilia-os/) (code) and [EOSVC](https://github.com/ersilia-os/eosvc) (data):

* Tracked by Git and linked to a GitHub repository: `scripts/`, `notebooks/`, `src/`, `tools/`, `docs/` and `assets/`.
* Tracked by `eosvc` and linked to a public or private S3 bucket: only `data/` and `output/`. The `access.json` file records whether they are public or private.
* Git LFS is not used. Files too large for Git (GitHub rejects files of 100 MiB or more) go in `data/` or `output/`, tracked by `eosvc`; keep files in `assets/` small.

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
