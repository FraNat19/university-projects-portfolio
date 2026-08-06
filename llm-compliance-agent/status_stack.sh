#!/bin/bash

set -euo pipefail

JOBID="${1:-}"

if [[ -z "${JOBID}" ]]; then echo "Usage: $0 JOBID"; exit 1; fi

squeue -j "${JOBID}" -o "%.18i %.9P %.20N %.8T %.10M %.20R"
