#!/bin/bash
#SBATCH --job-name=scrap_inail
#SBATCH --partition=low-wn
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/projects/CRISALIA/logs/slurm_inail_%j.out
#SBATCH --error=/projects/CRISALIA/logs/slurm_inail_%j.err
# ---------------------------------------------------------------------------
# Scraping COMPLETO catalogo INAIL (tutte le pagine: end_page 999 in hpc.yaml,
# lo scraper si ferma da solo a fine catalogo via max_empty_pages).
# Partizione low-wn: CPU+rete, nessuna GPU (verificato 16/07: i nodi di calcolo
# hanno accesso internet). Idempotente: salta gli ID gia' in data/metadata/.
# ---------------------------------------------------------------------------
set -uo pipefail
source /projects/CRISALIA/env.sh
cd /projects/CRISALIA/graphrag

echo "[job $SLURM_JOB_ID] nodo=$(hostname) inizio=$(date -Is)"
python -u src/scraping/inail/inail_scraper.py --config configs/hpc.yaml "$@"
rc=$?
echo "[job $SLURM_JOB_ID] fine=$(date -Is) exit=$rc"
exit $rc
