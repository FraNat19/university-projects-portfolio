#!/bin/bash
#SBATCH --job-name=recapH2
#SBATCH --partition=h100gpu
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=/projects/CRISALIA/logs/recapH2_%j.out
#SBATCH --error=/projects/CRISALIA/logs/recapH2_%j.err
set -uo pipefail
source /projects/CRISALIA/env.sh
cd /projects/CRISALIA/graphrag
echo "nodo=$(hostname) inizio=$(date -Is)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
python -u src/extraction/recaption_h2.py \
    --targets data/h2_targets.json \
    --max-new-tokens 512
echo "fine=$(date -Is)"
