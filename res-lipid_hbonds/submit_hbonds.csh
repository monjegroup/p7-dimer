#!/bin/sh
#SBATCH --ntasks-per-node=56
#SBATCH -N 1
#SBATCH --job-name=hbr3p1
#SBATCH -o hbr3p1.o%j
#SBATCH -t 24:00:00
#SBATCH --cluster=faculty
#SBATCH --partition=vmonje
#SBATCH --qos=vmonje
#SBATCH --mail-user=ocampbel@buffalo.edu
#SBATCH --mail-type=END

module load gcc/11.2.0 openmpi/4.1.1 vmd/1.9.4a57-CUDA-11.8.0

vmd -dispdev text -e calc_hbonds.tcl > calc_hbonds.out
