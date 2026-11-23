#!/usr/bin/env bash
# Explicit allocation: --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=32.
# Arguments: frozen benchmark_deepsomatic.py, prepared experiment directory.
set -eo pipefail
source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
module load StdEnv/2023 apptainer/1.3.5
srun --kill-on-bad-exit=1 --cpu-bind=cores python3 "$1" worker --out "$2"
