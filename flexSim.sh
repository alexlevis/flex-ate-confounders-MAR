#!/bin/bash
#SBATCH --job-name flexSim      # Set a name for your job. This is especially useful if you
#SBATCH --partition fasse     # Slurm partition to use
#SBATCH -c 1        # Number of tasks to run 3
#SBATCH --time 0-00:30       # Wall time limit in D-HH:MM
#SBATCH --mem 2000     # Memory limit for each tasks (in MB) # 1500
#SBATCH -o /n/home_fasse/awlevis/flex-sim/error/output_flexSim.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e /n/home_fasse/awlevis/flex-sim/error/errors_flexSim.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --array 1-5

# export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
module load R/4.2.2-fasrc01

Rscript '/n/home_fasse/awlevis/flex-sim/sim_temp.R'
