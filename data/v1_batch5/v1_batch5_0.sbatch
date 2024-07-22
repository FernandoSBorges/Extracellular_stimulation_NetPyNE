#!/bin/bash
#SBATCH --job-name=v1_batch5_0
#SBATCH -A TG-IBN140002
#SBATCH -t 23:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -o /expanse/lustre/projects/csd403/fborges/v1_batch5/v1_batch5_0.run
#SBATCH -e /expanse/lustre/projects/csd403/fborges/v1_batch5/v1_batch5_0.err
#SBATCH --mail-user=fernandodasilvaborges@gmail.com
#SBATCH --mail-type=end

#SBATCH --constraint="lustre"
#SBATCH --mem=1024G
#SBATCH --export=ALL
#SBATCH --partition=large-shared
source ~/.bashrc
cd /home/fborges/Extracellular_stimulation_NetPyNE/sim/
mpirun -n 128 nrniv -python -mpi init.py simConfig=/expanse/lustre/projects/csd403/fborges/v1_batch5/v1_batch5_0_cfg.json netParams=/expanse/lustre/projects/csd403/fborges/v1_batch5/v1_batch5_netParams.py
wait
    