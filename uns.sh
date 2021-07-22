#!/bin/bash -l

#SBATCH -o ./simu.out.%j
#SBATCH -e ./simu.err.%j

#SBATCH -D ./

#SBATCH -J simu
#SBATCH --partition=n0512
#SBATCH --nodes=37
#SBATCH --ntasks-per-node=40

#SBATCH --mem=80000
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=phe@ab.mpg.de
#SBATCH --time=24:00:00

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/cobra/u/penghe/local/lib"

module load openmpi/4
mpic++ -Wall -g uns.cpp -I/cobra/u/penghe/local/include/igraph -L/cobra/u/penghe/local/lib -ligraph -o uns
srun ./uns