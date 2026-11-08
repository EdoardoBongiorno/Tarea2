# TCGA-BRCA Differential Gene Analysis (PyDESeq2)

Edoardo Bongiorno
Lenguage: Python

---

# Project Overview
This project implements a differential gene expression analysis pipeline on the TCGA-BRCA dataset (Tumor vs. Normal tissue). The workflow translates the original R/DESeq2 pipeline into Python using PyDESeq2, maintaining statistical equivalence, identical filtering criteria, and biological coherence as specified in the course guidelines.

---

# Repository Structure & GitHub Link
GitHub repository: https://github.com/EdoardoBongiorno/Tarea2.git
Directory Layout:
 `src/`: Contains the main Python execution script (`01_run_deseq2.py`).
 `data/raw/`: Contains input files (`counts.csv.gz`, `metadata.csv.gz`). *Note: Excluded from GitHub.*
 `results/`: Contains output TSV tables and PNG figures. *Note: Excluded from GitHub.*
 `logs/`: Contains execution summaries and package version logs.

---

# Instructions to Reproduce the Analysis

Create and Activate the Conda Environment
The project environment is defined in the `environment.yml` file. Run the following command in your terminal to build it:
```bash
conda env create -f environment.yml
conda activate pydeseq2