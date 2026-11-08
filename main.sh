#!/bin/bash
cd "$(dirname "$0")"
eval "$(conda shell.bash hook)"
conda activate pydeseq2
python src/01_run_deseq2.py