#!/bin/bash

set -euo pipefail

JOBID=$(sbatch --parsable /mnt/beegfs/proj/dss.dmaia/thesis_graph_rag/slurm_stack_qdrant_neo4j.sbatch)

echo "JOBID=${JOBID}"

squeue -j "${JOBID}" -o "%.18i %.20N %.8T %.10M"
