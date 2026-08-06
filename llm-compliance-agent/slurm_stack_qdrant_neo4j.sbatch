#!/bin/bash

#SBATCH --job-name=graphrag_stack

#SBATCH --partition=low-prio

#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --gres=gpu:1

#SBATCH --cpus-per-task=8

#SBATCH --mem=32G

#SBATCH --time=24:00:00

#SBATCH --output=slurm_stack_%j.log

#SBATCH --error=slurm_stack_%j.log



set -euo pipefail



BASE="/mnt/beegfs/proj/dss.dmaia/thesis_graph_rag"



QDRANT_SIF="${BASE}/qdrant.sif"

NEO4J_SIF="${BASE}/neo4j.sif"



LOGDIR="${BASE}/logs"

QDRANT_STORAGE="${BASE}/qdrant_storage_job_${SLURM_JOB_ID}"



NEO4J_DATA="${BASE}/neo4j_data"

NEO4J_LOGS="${BASE}/neo4j_logs"

NEO4J_CONF="${BASE}/neo4j_conf"

NEO4J_PLUGINS="${BASE}/neo4j_plugins"

NEO4J_IMPORT="${BASE}/neo4j_import"



echo "HOSTNAME=$(hostname)"

echo "JOBID=${SLURM_JOB_ID}"

echo "BASE=${BASE}"

echo "QDRANT_SIF=${QDRANT_SIF}"

echo "NEO4J_SIF=${NEO4J_SIF}"



ls -lh "${QDRANT_SIF}" "${NEO4J_SIF}"



module load apptainer

echo "APPTAINER=$(command -v apptainer)"

apptainer --version



mkdir -p "${LOGDIR}" "${QDRANT_STORAGE}" "${NEO4J_DATA}" "${NEO4J_LOGS}" "${NEO4J_CONF}" "${NEO4J_PLUGINS}" "${NEO4J_IMPORT}"

chmod 777 "${QDRANT_STORAGE}" || true



export NEO4J_AUTH="neo4j/${NEO4J_PASSWORD:?set NEO4J_PASSWORD before submitting}"

export NEO4J_server_default__listen__address="0.0.0.0"

export NEO4J_server_http_listen__address="0.0.0.0:7474"

export NEO4J_server_https_enabled="false"

export NEO4J_server_bolt_listen__address="0.0.0.0:7687"



echo "STARTING_QDRANT"

apptainer run --writable-tmpfs --bind "${QDRANT_STORAGE}:/qdrant/storage" "${QDRANT_SIF}" > "${LOGDIR}/qdrant_${SLURM_JOB_ID}.log" 2>&1 &

QDRANT_PID=$!

echo "QDRANT_PID=${QDRANT_PID}"



echo "STARTING_NEO4J"

apptainer run --writable-tmpfs --bind "${NEO4J_DATA}:/data" --bind "${NEO4J_LOGS}:/logs" --bind "${NEO4J_CONF}:/conf" --bind "${NEO4J_PLUGINS}:/plugins" --bind "${NEO4J_IMPORT}:/import" "${NEO4J_SIF}" > "${NEO4J_LOGS}/neo4j_${SLURM_JOB_ID}.log" 2>&1 &

NEO4J_PID=$!

echo "NEO4J_PID=${NEO4J_PID}"



# Wait readiness

for i in $(seq 1 120); do

  Q_OK=0; N_OK=0

  curl -sS --max-time 1 http://localhost:6333/collections >/dev/null 2>&1 && Q_OK=1

  (echo > /dev/tcp/127.0.0.1/7687) >/dev/null 2>&1 && N_OK=1

  if [[ "${Q_OK}" -eq 1 && "${N_OK}" -eq 1 ]]; then

    echo "STACK READY"

    break

  fi

  sleep 1

done



# Hard fail with logs if not ready

curl -sS --max-time 1 http://localhost:6333/collections >/dev/null 2>&1 || { echo "Qdrant not ready"; tail -n 120 "${LOGDIR}/qdrant_${SLURM_JOB_ID}.log" || true; exit 1; }

(echo > /dev/tcp/127.0.0.1/7687) >/dev/null 2>&1 || { echo "Neo4j not ready"; tail -n 200 "${NEO4J_LOGS}/neo4j_${SLURM_JOB_ID}.log" || true; exit 1; }



echo "Keep-alive..."

wait

