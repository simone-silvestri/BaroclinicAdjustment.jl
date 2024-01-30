#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --mem=200GB
#SBATCH --time=120:00:00
#SBATCH -p sched_mit_raffaele_gpu

source setup_engaging.sh

export CASE="Upwind"

$JULIA --project --check-bounds=no -e 'using Pkg; Pkg.instantiate()'
$JULIA --project --check-bounds=no run_restoring.jl
