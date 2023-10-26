#!/bin/bash
#SBATCH --job-name flexSim_par      # Set a name for your job. This is especially useful if you
#SBATCH --partition fasse     # Slurm partition to use
#SBATCH -c 1        # Number of tasks to run 1
#SBATCH --time 0-00:07       # Wall time limit in D-HH:MM
#SBATCH --mem 800     # Memory limit for each tasks (in MB) # 800
#SBATCH -o /n/home_fasse/awlevis/flex-sim/error/output_flexSim_par.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e /n/home_fasse/awlevis/flex-sim/error/errors_flexSim_par.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --array 1-2500

# export R_LIBS_USER=$HOME/apps/R_3.5.1:$R_LIBS_USER
module load R/4.2.2-fasrc01

Rscript '/n/home_fasse/awlevis/flex-sim/sim_temp.R'
