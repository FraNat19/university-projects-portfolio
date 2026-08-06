#!/bin/bash
#SBATCH --job-name=inail_enrichment
#SBATCH --output=/mnt/beegfs/proj/dss.dmaia/INAILProj/logs/slurm_%j.out
#SBATCH --error=/mnt/beegfs/proj/dss.dmaia/INAILProj/logs/slurm_%j.err
#SBATCH --partition=defq-noprio
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# ============================================================================
# SLURM Job Script per INAIL Document Enrichment
# ============================================================================

source /mnt/beegfs/home/fnatali/docling/inail_env/bin/activate

echo "========================================="
echo "INAIL Enrichment Job Started"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "========================================="

# Carica moduli necessari
# Nessun module: uso direttamente python3 del sistema
# Se in futuro crei un venv, attivalo così:
# source /mnt/beegfs/home/fnatali/docling/inail_env/bin/activate

# Verifica GPU disponibile
nvidia-smi

# Directory di lavoro
cd /mnt/beegfs/proj/dss.dmaia/INAILProj

# Esegui script Python
echo ""

echo "Avvio elaborazione documenti..."

python3 /mnt/beegfs/home/fnatali/docling/tf_idf_embedgemma_hpc.py

# Capture exit code
EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job Completed"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "========================================="

exit $EXIT_CODE
