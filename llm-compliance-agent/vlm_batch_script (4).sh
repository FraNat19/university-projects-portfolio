#!/bin/bash

#SBATCH --job-name=vlm_images
#SBATCH --output=logs/vlm_%j.out
#SBATCH --error=logs/vlm_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=defq-noprio     
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --gres=gpu:a100:1

# ===== CONFIGURAZIONE =====

# Modifica questi path con i tuoi

IMAGES_DIR="/mnt/beegfs/home/fnatali/docling/data/images"

JSON_DIR="/mnt/beegfs/home/fnatali/docling/data/json"

VENV_PATH="/mnt/beegfs/home/fnatali/docling/venv"  # Se hai un venv

SCRIPT_PATH="./apply_vlm_to_existing_images.py"



# Parametri VLM (modificabili da command line)

BATCH_SIZE=${BATCH_SIZE:-8}

IMAGE_SIZE=${IMAGE_SIZE:-512}

USE_8BIT=${USE_8BIT:-false}

FORCE=${FORCE:-false}



# ===== SETUP =====

echo "========================================="

echo "VLM APPLICATION JOB"

echo "========================================="

echo "Job ID: $SLURM_JOB_ID"

echo "Node: $SLURM_NODELIST"

echo "Start time: $(date)"

echo "GPU: $CUDA_VISIBLE_DEVICES"

echo "========================================="

echo ""


# Nella sezione moduli:
module purge
module load gcc
module load python/3.12.1
module load cuda/11.8.0

echo "Moduli caricati:"

module list

echo ""



# Attiva ambiente virtuale se esiste

if [ -d "$VENV_PATH" ]; then

    echo "Attivazione ambiente virtuale..."

    source "${VENV_PATH}/bin/activate"

    echo "Python: $(which python)"

    echo "Python version: $(python --version)"

else

    echo "[WARNING] Ambiente virtuale non trovato, uso Python di sistema"

fi

echo ""



# Verifica GPU

echo "Verifica GPU..."

nvidia-smi --query-gpu=name,memory.total,memory.free --format=csv

echo ""



# ===== PARAMETRI =====

echo "========================================="

echo "PARAMETRI VLM"

echo "========================================="

echo "Images directory: $IMAGES_DIR"

echo "JSON directory: $JSON_DIR"

echo "Batch size: $BATCH_SIZE"

echo "Image size: ${IMAGE_SIZE}px"

echo "8-bit quantization: $USE_8BIT"

echo "Force reprocess: $FORCE"

echo "========================================="

echo ""



# Verifica directory

if [ ! -d "$IMAGES_DIR" ]; then

    echo "[ERROR] Images directory non esiste: $IMAGES_DIR"

    exit 1

fi



# Conta immagini

TOTAL_PNG=$(find "$IMAGES_DIR" -name "*.png" | wc -l)

TOTAL_TXT=$(find "$IMAGES_DIR" -name "*.txt" | wc -l)

TO_PROCESS=$((TOTAL_PNG - TOTAL_TXT))



echo "Statistiche pre-processing:"

echo "  - Immagini PNG totali: $TOTAL_PNG"

echo "  - Descrizioni esistenti: $TOTAL_TXT"

echo "  - Da processare: $TO_PROCESS"

echo ""



if [ "$TO_PROCESS" -le 0 ] && [ "$FORCE" != "true" ]; then

    echo "[INFO] Tutte le immagini hanno già descrizioni!"

    echo "[INFO] Usa FORCE=true per riprocessare"

    exit 0

fi



# ===== ESECUZIONE =====

START_TIME=$(date +%s)



# Costruisci comando

CMD="python $SCRIPT_PATH \

    --images-dir $IMAGES_DIR \

    --json-dir $JSON_DIR \

    --batch-size $BATCH_SIZE \

    --image-size $IMAGE_SIZE"



if [ "$USE_8BIT" = "true" ]; then

    CMD="$CMD --use-8bit"

fi



if [ "$FORCE" = "true" ]; then

    CMD="$CMD --force"

fi



echo "Comando: $CMD"

echo ""

echo "Inizio processing..."

echo "========================================="

echo ""



# Esegui con auto-conferma

echo "y" | eval $CMD



EXIT_CODE=$?



END_TIME=$(date +%s)

DURATION=$((END_TIME - START_TIME))

HOURS=$((DURATION / 3600))

MINUTES=$(((DURATION % 3600) / 60))

SECONDS=$((DURATION % 60))



# ===== REPORT =====

echo ""

echo "========================================="

echo "JOB COMPLETATO"

echo "========================================="

echo "Exit code: $EXIT_CODE"

echo "End time: $(date)"

echo "Durata: ${HOURS}h ${MINUTES}m ${SECONDS}s"

echo ""



# Statistiche post-processing

TOTAL_TXT_AFTER=$(find "$IMAGES_DIR" -name "*.txt" | wc -l)

CREATED=$((TOTAL_TXT_AFTER - TOTAL_TXT))



echo "Statistiche post-processing:"

echo "  - Descrizioni totali: $TOTAL_TXT_AFTER"

echo "  - Nuove create: $CREATED"

echo ""

echo "========================================="



exit $EXIT_CODE
