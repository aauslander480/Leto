#!/bin/bash
#SBATCH --job-name=leto-qgis-pilot
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=leto-qgis.%j.out
#SBATCH --error=leto-qgis.%j.err

set -euo pipefail
# Export these before sbatch. Find your exact module with: ml spider QGIS
: "${LETO_QGIS_MODULE:?Set LETO_QGIS_MODULE to the installed QGIS module/version}"
: "${LETO_DATA_ROOT:?Set LETO_DATA_ROOT to your transferred Maine data}"
: "${LETO_RUNS_ROOT:?Set LETO_RUNS_ROOT to your scratch output directory}"
module load "$LETO_QGIS_MODULE"
export QT_QPA_PLATFORM=offscreen
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
cd "$SLURM_SUBMIT_DIR"
# Set QGIS_PYTHON if your module uses a different Python executable.
"${QGIS_PYTHON:-python3}" run_leto.py --check-runtime
"${QGIS_PYTHON:-python3}" run_leto.py configs/maine_sapelo2.py
