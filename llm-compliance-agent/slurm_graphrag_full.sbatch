#!/bin/bash
#SBATCH --job-name=graphrag_full
#SBATCH --partition=low-prio
#SBATCH --nodes=1
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/logs/graphrag_full_%j.log

set -euo pipefail

echo "========================================"
echo "🚀 GRAPHRAG FULL STACK - ALL IN ONE"
echo "========================================"
echo "🕐 Started: $(date)"
echo "📍 Node: $(hostname)"
echo "📋 Job ID: ${SLURM_JOB_ID}"
echo ""

BASE="/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag"
cd "$BASE"

# ================== GPU INFO ==================
echo "🎮 GPU INFO:"
nvidia-smi --query-gpu=name,memory.total --format=csv
echo ""

# ================== LOAD APPTAINER ==================
module load apptainer
echo "✅ Apptainer loaded"

# ================== DIRECTORIES ==================
QDRANT_SIF="${BASE}/qdrant.sif"
NEO4J_SIF="${BASE}/neo4j.sif"
LOGDIR="${BASE}/logs"
QDRANT_STORAGE="${BASE}/qdrant_storage"
NEO4J_DATA="${BASE}/neo4j_data"
NEO4J_LOGS="${BASE}/neo4j_logs"
NEO4J_CONF="${BASE}/neo4j_conf"
NEO4J_PLUGINS="${BASE}/neo4j_plugins"
NEO4J_IMPORT="${BASE}/neo4j_import"

mkdir -p "${LOGDIR}" "${QDRANT_STORAGE}" "${NEO4J_DATA}" "${NEO4J_LOGS}" "${NEO4J_CONF}" "${NEO4J_PLUGINS}" "${NEO4J_IMPORT}"

# ================== START QDRANT ==================
echo ""
echo "📦 Starting Qdrant..."
apptainer run --writable-tmpfs \
    --bind "${QDRANT_STORAGE}:/qdrant/storage" \
    "${QDRANT_SIF}" > "${LOGDIR}/qdrant_${SLURM_JOB_ID}.log" 2>&1 &
QDRANT_PID=$!
echo "   PID: $QDRANT_PID"

# ================== START NEO4J ==================
echo "📊 Starting Neo4j..."
export NEO4J_AUTH="neo4j/${NEO4J_PASSWORD:?set NEO4J_PASSWORD before submitting}"
export NEO4J_server_default__listen__address="0.0.0.0"
export NEO4J_server_http_listen__address="0.0.0.0:7474"
export NEO4J_server_bolt_listen__address="0.0.0.0:7687"

apptainer run --writable-tmpfs \
    --bind "${NEO4J_DATA}:/data" \
    --bind "${NEO4J_LOGS}:/logs" \
    --bind "${NEO4J_CONF}:/conf" \
    --bind "${NEO4J_PLUGINS}:/plugins" \
    --bind "${NEO4J_IMPORT}:/import" \
    "${NEO4J_SIF}" > "${NEO4J_LOGS}/neo4j_${SLURM_JOB_ID}.log" 2>&1 &
NEO4J_PID=$!
echo "   PID: $NEO4J_PID"

# ================== START OLLAMA ==================
echo "🤖 Starting Ollama with Qwen 72B..."
export OLLAMA_MODELS="${BASE}/ollama_models_72b"
export OLLAMA_HOST="127.0.0.1:11434"

./ollama serve > "${LOGDIR}/ollama_${SLURM_JOB_ID}.log" 2>&1 &
OLLAMA_PID=$!
echo "   PID: $OLLAMA_PID"

# ================== WAIT FOR SERVICES ==================
echo ""
echo "⏳ Waiting for services to start (90 seconds max)..."

for i in $(seq 1 90); do
    QDRANT_OK=0
    NEO4J_OK=0
    OLLAMA_OK=0
    
    # Check Qdrant
    curl -sS --max-time 2 http://localhost:6333/collections > /dev/null 2>&1 && QDRANT_OK=1
    
    # Check Neo4j
    (echo > /dev/tcp/127.0.0.1/7687) > /dev/null 2>&1 && NEO4J_OK=1
    
    # Check Ollama
    curl -sS --max-time 2 http://localhost:11434/api/tags > /dev/null 2>&1 && OLLAMA_OK=1
    
    echo -n "   [$i/90] Qdrant:$QDRANT_OK Neo4j:$NEO4J_OK Ollama:$OLLAMA_OK"
    
    if [[ "$QDRANT_OK" -eq 1 && "$NEO4J_OK" -eq 1 && "$OLLAMA_OK" -eq 1 ]]; then
        echo " ✅ ALL READY!"
        break
    fi
    echo ""
    sleep 1
done

# ================== VERIFY ==================
echo ""
echo "🔍 Verifying services..."

# Qdrant
if curl -sS http://localhost:6333/collections > /dev/null 2>&1; then
    COLLECTIONS=$(curl -s http://localhost:6333/collections | python3 -c "import sys,json; print(len(json.load(sys.stdin).get('result',{}).get('collections',[])))" 2>/dev/null || echo "?")
    echo "   ✅ Qdrant: $COLLECTIONS collections"
else
    echo "   ❌ Qdrant NOT responding"
    tail -20 "${LOGDIR}/qdrant_${SLURM_JOB_ID}.log"
    exit 1
fi

# Neo4j
if (echo > /dev/tcp/127.0.0.1/7687) > /dev/null 2>&1; then
    echo "   ✅ Neo4j: Bolt ready on :7687"
else
    echo "   ❌ Neo4j NOT responding"
    tail -20 "${NEO4J_LOGS}/neo4j_${SLURM_JOB_ID}.log"
    exit 1
fi

# Ollama
if curl -s http://localhost:11434/api/tags > /dev/null 2>&1; then
    MODELS=$(./ollama list 2>/dev/null | tail -n +2 | wc -l)
    echo "   ✅ Ollama: $MODELS models available"
    ./ollama list
else
    echo "   ❌ Ollama NOT responding"
    tail -20 "${LOGDIR}/ollama_${SLURM_JOB_ID}.log"
    exit 1
fi

# ================== ENVIRONMENT FOR PYTHON ==================
echo ""
echo "📝 Setting environment variables..."

# Per Python - usa porte LOCALI (non tunnel!)
export QDRANT_HOST="127.0.0.1"
export QDRANT_PORT="6333"
export NEO4J_URI="bolt://127.0.0.1:7687"
export NEO4J_USER="neo4j"
export NEO4J_PASSWORD="${NEO4J_PASSWORD:?set NEO4J_PASSWORD before submitting}"
export OLLAMA_URL="http://127.0.0.1:11434"

echo "   QDRANT_PORT=$QDRANT_PORT"
echo "   NEO4J_URI=$NEO4J_URI"
echo "   OLLAMA_URL=$OLLAMA_URL"

# ================== READY ==================
echo ""
echo "========================================"
echo "✅ GRAPHRAG FULL STACK READY!"
echo "========================================"
echo ""
echo "📊 Services running on $(hostname):"
echo "   - Qdrant:  http://localhost:6333"
echo "   - Neo4j:   bolt://localhost:7687 (HTTP: :7474)"
echo "   - Ollama:  http://localhost:11434"
echo ""
echo "🤖 Model: Qwen 2.5 72B"
echo "🎮 GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo ""
echo "📂 Logs: ${LOGDIR}/"
echo ""
echo "🚀 To run the RAG system, use srun --overlap from another terminal:"
echo ""
echo "   srun --overlap --jobid=${SLURM_JOB_ID} python3 graphrag_rag_complete.py -i"
echo ""
echo "   Or modify the script below to run automatically."
echo ""
echo "========================================"

# ================== OPTIONAL: RUN RAG AUTOMATICALLY ==================
# Uncomment these lines to run the RAG system automatically:

# echo ""
# echo "🎬 Running GraphRAG interactive mode..."
# source venv/bin/activate
# python3 graphrag_rag_complete.py --interactive \
#     --qdrant-port 6333 \
#     --neo4j-port 7687 \
#     --ollama-host http://127.0.0.1:11434

# ================== KEEP ALIVE ==================
echo ""
echo "⏳ Keeping services alive..."
echo "   Press Ctrl+C to stop (or job will end at time limit)"
echo ""

# Wait for all background processes
wait
