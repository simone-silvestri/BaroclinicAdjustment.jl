# Upload modules
source /etc/profile.d/modules.sh
module load cuda/11.3
# module load openmpi/4.1.2
# module load gcc/12.2.0-x86_64
# module load openmpi/4.1.4-pmi-cuda-ucx-x86_64
# module load nvhpc/2023_233/nvhpc/23.3

export JULIA_CUDA_USE_BINARYBUILD=false
export JULIA_CUDA_MEMORY_POOL=none
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_NVTX_CALLBACKS=gc
export JULIA_DEPOT_PATH="/orcd/nese/raffaele/001/ssilvest/depot"
export JULIA="/nfs/cnhlab001/ssilvest/julia-1.9.2/bin/julia"

