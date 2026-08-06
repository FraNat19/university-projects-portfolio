#!/bin/bash
#SBATCH --job-name=scrap_osha
#SBATCH --partition=low-wn
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/projects/CRISALIA/logs/slurm_osha_%j.out
#SBATCH --error=/projects/CRISALIA/logs/slurm_osha_%j.err
# ---------------------------------------------------------------------------
# Scraping OSHA prime 100 pagine (start_page 0 -> end_page 99 in hpc.yaml,
# decisione utente 2026-07-05). Lingua PDF: EN > IT.
# Partizione low-wn: CPU+rete, nessuna GPU (verificato 16/07: i nodi di calcolo
# hanno accesso internet). Idempotente: salta gli ID gia' in data/metadata/.
# ---------------------------------------------------------------------------
set -uo pipefail
source /projects/CRISALIA/env.sh
cd /projects/CRISALIA/graphrag

echo "[job $SLURM_JOB_ID] nodo=$(hostname) inizio=$(date -Is)"
python -u src/scraping/osha/osha_scraper.py --config configs/hpc.yaml "$@"
rc=$?
echo "[job $SLURM_JOB_ID] fine=$(date -Is) exit=$rc"
exit $rc
