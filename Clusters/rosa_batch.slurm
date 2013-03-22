#!/bin/bash
#SBATCH --job-name="FANCYNAME"
#SBATCH --ntasks=NCPU
#SBATCH --time=9:00:00

aprun -n NCPU ./MY_SIMULATION_BINARY
