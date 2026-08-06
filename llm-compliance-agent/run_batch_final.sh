#!/bin/bash

#SBATCH --job-name=inail_batch

#SBATCH --output=logs/batch_%j.out

#SBATCH --error=logs/batch_%j.err

#SBATCH --time=72:00:00

#SBATCH --partition=defq-noprio

#SBATCH --nodes=1

#SBATCH --cpus-per-task=8

#SBATCH --mem=64GB



echo "========================================" 

echo "INAIL Batch Scraping - Job $SLURM_JOB_ID"

echo "Node: $SLURM_NODELIST"

echo "Start: $(date)"

echo "========================================"



# Export variabili per headless mode

export DISPLAY=:99

export CHROME_BIN=/usr/bin/google-chrome

export MOZ_HEADLESS=1



# Setup

cd /mnt/beegfs/home/fnatali/docling

module load gcc python

source venv/bin/activate



echo "Python: $(which python)"

echo "Starting scraper..."



# Lancia con stdbuf per vedere output in tempo reale

stdbuf -oL -eL python -u inail_hpc_main.py \

  --start-page 1 \

  --end-page 87 \

  --max-docs 864 \

  --save-pdfs \

  --no-vlm \

  --output-dir /mnt/beegfs/home/fnatali/docling/data



EXIT_CODE=$?



echo "========================================"

echo "Exit code: $EXIT_CODE"

echo "End: $(date)"

echo "========================================"



exit $EXIT_CODE

