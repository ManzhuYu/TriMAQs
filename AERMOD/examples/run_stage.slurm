#!/bin/bash
#SBATCH --job-name=TriMAQs_AERMOD
#SBATCH --output=TriMAQs_%j.out
#SBATCH --error=TriMAQs_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=160G

# Set the account, partition and environment for your HPC system.
set -eo pipefail
# module load anaconda
# eval "$(conda shell.bash hook)"
# conda activate YOUR_ENVIRONMENT

CODE_ROOT="/path/to/AERMOD"
python -u "$CODE_ROOT/run.py" run_production --config "$CODE_ROOT/config/settings.json"
